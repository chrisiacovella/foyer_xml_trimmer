"""
This module will generate a Foyer XML file that contains
only the subset of parameters that appear in an atom typed molecule.

To ensure that the trimmed XML and original XML file can be directly compared
(i.e., to ensure no information is lost or reordered in topological parameters)
this function will take as input both a typed Parmed structure object
and a Foyer XML file.

"""

import mbuild as mb
import parmed as pmd
import foyer

from foyer import Forcefield

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element,tostring, ElementTree
from pkg_resources import resource_filename

# returns an XML Element for a user defined tag (str)
# with attributes that are defined in temp_dict (dictionary)
def _dict_to_xml(tag, temp_dict):
 
    elem = Element(tag)
    elem.attrib = temp_dict
         
    return elem

# topological parameters sets can be defined for a collection of atoms by class,
# e.g., class1, class2, class3, class4 (for a dihedral)
# or for a collection of atoms as defined by atomtype
# e.g., type1, type2, type3, type4, or a mixture therein
# e.g., type1, class2, class3, type4
# This function will look at the attributes in attrib_dict
# to identify what has been defined in the XML so we can do the appropriate matching.
# This defines a weight to each schema, where least specific have the largest weight
# i.e., each time we encounter class we update the weight by 1
# and most specific have the lowest weight value, i.e., they add 0 to the running weight.
def _identify_schema(attrib_dict, nparams=2):
    schema = []
    weight = 0
    for i in range(1, nparams+1):
        if f'class{i}' in attrib_dict:
            schema.append(f'class{i}')
            weight +=1
        else:
            schema.append(f'type{i}')
    
    return schema, weight

# simple helper function that will return the class type if class appears in the schema
# otherwise the atom_type is returned
def _switch_class_type(schema, atom_class, atom_type):
    if 'class' in schema:
        return atom_class
    else:
        return atom_type


# Function that will search the XML data for the equivalent topological collection and write the parameters to
# the new xml file output.
def _topology_match(atom_type_dict, typed_topo, xml_root, blank_root, topo_type, n_params):

        #define a dictionary where the val corresponds to the right section in the output xml
        xml_section = {
                        'Bond': 'HarmonicBondForce',
                        'Angle':'HarmonicAngleForce',
                        'Proper': 'RBTorsionForce',
                        'Improper': 'PeriodicTorsionForce'
                    }
        #define a dictionary that will point us to the correct name in parmed data structure
        parmed_section = {
                        'Bond': 'bonds',
                        'Angle':'angles',
                        'Proper': 'rb_torsions',
                        'Improper': 'impropers'
                    }
        #definition of equivalent sequences for parameter matching
        #bond/angle/proper all follows simple forward/reverse schemes,
        #impropers are defined about a center atom, which appears first in the list,
        #and thus any order of the other 3 atoms is equivalent
        sequences = {
                        'Bond': [[0,1], [1,0]],
                        'Angle': [[0,1,2],[2,1,0]],
                        'Proper': [[0,1,2,3], [3,2,1,0]],
                        'Improper': [ [0,1,2,3], [0,1,3,2], [0,2,1,3], [0,2,3,1], [0,3,1,2], [0,3,2,1]],
        
                    }
        topo_by_type = []

        # rather than using a lot of if/else statements, just use eval to get us
        # into the appropriate part of the parmed data structure.
        topo_temp = eval(f'typed_topo.{parmed_section[topo_type]}')
        for topo_element in topo_temp:
            
            # use eval again rather than lots of if/else statements
            topo_af = []
            for i in range(0,n_params):
                topo_af.append(eval(f'topo_element.atom{i+1}.type'))
            
            # We don't want to include duplicates, including equivalent sequences.
            # We will check if this collection of atom_types or equivalent sequences
            # already have been added to in the topo_type array
            unique = True
            for seq in sequences[topo_type]:
                if [topo_af[i] for i in seq] in topo_by_type:
                    unique = False
                    break
            if unique:
                topo_by_type.append(topo_af)
        
        # entries can be defined using class, atom_type, or a combination of both
        # this code will examine each entry in the xml file and identifies the schema
        # that way we make the right comparison later
        topo_entry_temp = []
        for topo in xml_root.iter(topo_type):
            schema, weight = _identify_schema(topo.attrib, n_params)
            topo_dict = {'schema': schema, 'weight': weight, 'attrib': topo.attrib}
            topo_entry_temp.append(topo_dict)
        
        # sort the list such that the most specific entries will be checked first.
        # Weights are assigned such that the least specific have the highest value
        # so we can sort from lowest to highest.
        topo_entries = sorted(topo_entry_temp, key=lambda d: d['weight'])
        
        unique_collection = []
        
        for topo_element in topo_by_type:
            not_matched = True
            for topo_entry in topo_entries:
                schema = topo_entry['schema']
                topo = topo_entry['attrib']

                for sequence in sequences[topo_type]:
                    match = []
                    collection = []
                    
                    if not_matched:
                        seq_test = [topo_element[i] for i in sequence]
                        for i,seq in enumerate(seq_test):
                            match.append(False)
                            if topo[schema[i]] == _switch_class_type(schema[i], atom_type_dict[seq] ,seq):
                                match[i] = True
                            collection.append(topo[schema[i]])
                                
                        if all(match):
                            if collection not in unique_collection:
                                unique_collection.append(collection)
                                elem = _dict_to_xml(topo_type, topo)
                                for child in blank_root:
                                    if child.tag == xml_section[topo_type]:
                                        child.append(elem)
                            not_matched = False
                            
def forcefield_trim(typed_molecule, input_xml, output_xml):

    # Read in an empty Foyer xml file with all the expected elements
    path_to_blank_xml  = resource_filename("foyer_xml_trimmer", "data/blank.xml")
    blank_tree = ET.parse(path_to_blank_xml)
    blank_root = blank_tree.getroot()
    
    # Read in the untrimmed foyer XML file
    xml_tree = ET.parse(input_xml)
    xml_root = xml_tree.getroot()
       
    # we need to check to ensure that the typed_molecule
    # is in fact a parmed structure object, otherwise
    # raise an error
    if not isinstance(typed_molecule, pmd.Structure):
        raise ValueError(
                    f"Expected Parmed Structure object, instead typed_molecule "
                    f"found to be {type(typed_molecule)}."
                )
    else:

        # Initialize a dict containing the atom_types
        # as the keys.  The associated value will be set
        # to an empty string initially. This value
        # will be later assigned info from the xml file
        atom_type_dict = {}

        for atom in typed_molecule.atoms:
            if atom.atom_type.name not in atom_type_dict:
                atom_type_dict[atom.atom_type.name] = ''
        
        #loop over each atom_type to set class and find if there are any overrides we need to define
        iterate = True
        while iterate:
            atom_type_overrides_dict = {}
            for atom_type in atom_type_dict:
                #loop over the Type definitions in the xml file
                for atom in xml_root.iter('Type'):
                    #find the matching atom_type
                    if atom.attrib['name'] == atom_type:
                        #associated the atom-class with atom-type
                        atom_type_dict[atom_type] = atom.attrib['class']
                        if 'overrides' in atom.attrib:
                            overrides = atom.attrib['overrides']
                            overrides = overrides.split(',')
                            for override in overrides:
                                if override not in atom_type_dict:
                                    atom_type_overrides_dict[override] = ''
            atom_type_dict.update(atom_type_overrides_dict)
        
            if len(atom_type_overrides_dict) > 0:
                for key in atom_type_overrides_dict:
                    print('Note: atom type: ', key, ' is referenced in an overrides statement, but does not appear in the system.')
            else:
                iterate = False
                
        # first we will loop over the AtomTypes
        for atom in xml_root.iter('Type'):
            if atom.attrib['name'] in atom_type_dict:
                atom_type_dict[atom.attrib['name']] = atom.attrib['class']
                elem = _dict_to_xml('Type', atom.attrib)
                
                #populate the blank XML file
                for child in blank_root:
                    if child.tag == 'AtomTypes':
                        child.append(elem)
            
        # next, loop over the NonBondedForce section
        for atom in xml_root.iter('Atom'):
            if atom.attrib['type'] in atom_type_dict:
                elem = _dict_to_xml('Atom', atom.attrib)
                
                #populate the blank XML file
                for child in blank_root:
                    if child.tag == 'NonbondedForce':
                        child.append(elem)
        
        _topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Bond', n_params=2)
        _topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Angle', n_params=3)
        _topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Proper', n_params=4)
        _topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Improper', n_params=4)
 
        ET.indent(blank_root, space='\t', level=0)
        blank_tree.write(output_xml, encoding="utf-8", xml_declaration=True)


def forcefield_score(typed_molecule, input_xml):

    
    # Read in the untrimmed foyer XML file
    xml_tree = ET.parse(input_xml)
    xml_root = xml_tree.getroot()
       
    # we need to check to ensure that the typed_molecule
    # is in fact a parmed structure object, otherwise
    # raise an error
    if not isinstance(typed_molecule, pmd.Structure):
        raise ValueError(
                    f"Expected Parmed Structure object, instead typed_molecule "
                    f"found to be {type(typed_molecule)}."
                )
    else:

        # Initialize a dict containing the atom_types
        # as the keys.  The associated value will be set
        # to an empty string initially. This value
        # will be later assigned info from the xml file
        atom_type_dict = {}

        for atom in typed_molecule.atoms:
            if atom.atom_type.name not in atom_type_dict:
                atom_type_dict[atom.atom_type.name] = ''
        
        #loop over each atom_type to set class and find if there are any overrides we need to define
        iterate = True
        while iterate:
            atom_type_overrides_dict = {}
            for atom_type in atom_type_dict:
                #loop over the Type definitions in the xml file
                for atom in xml_root.iter('Type'):
                    #find the matching atom_type
                    if atom.attrib['name'] == atom_type:
                        #associated the atom-class with atom-type
                        atom_type_dict[atom_type] = atom.attrib['class']
                        if 'overrides' in atom.attrib:
                            overrides = atom.attrib['overrides']
                            overrides = overrides.split(',')
                            for override in overrides:
                                if override not in atom_type_dict:
                                    atom_type_overrides_dict[override] = ''
            atom_type_dict.update(atom_type_overrides_dict)
        
            if len(atom_type_overrides_dict) > 0:
                for key in atom_type_overrides_dict:
                    print('Note: atom type: ', key, ' is referenced in an overrides statement, but does not appear in the system.')
            else:
                iterate = False
        xml_atom_type_dict = {}
        xml_atom_type_dict_score = {}
        # first we will loop over the AtomTypes
        for atom in xml_root.iter('Type'):
            xml_atom_type_dict[atom.attrib['name']] = atom.attrib['class']
            xml_atom_type_dict_score[atom.attrib['name']] = 0
            
        for atom in xml_atom_type_dict:
            if atom in atom_type_dict:
                xml_atom_type_dict_score
                
        
        #_topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Bond', n_params=2)
        #_topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Angle', n_params=3)
        #_topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Proper', n_params=4)
        #_topology_match(atom_type_dict=atom_type_dict, typed_topo=typed_molecule, xml_root=xml_root, blank_root=blank_root, topo_type='Improper', n_params=4)
 




