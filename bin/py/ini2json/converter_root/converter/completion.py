# -*- coding: utf-8 -*-
'''
Created on 15.4.2012

@author: Tomáš Košek
'''
import materials as mtr #@UnresolvedImport

import json
import re

import sys, os 

def complete (adr, res, warn, comm, dbg, alt, mat):
    """
    Starts the second convertion of specified file.
    Converted data prints in specified result file.

    :param lat: adr
    :type lat: string
    :param lat: res
    :type lat: string
    :param lat: warn
    :type lat: bool
    :param lat: comm
    :type lat: bool
    :param lat: debug
    :type lat: bool
    :param lat: alt
    :type lat: bool
    :param lat: mat
    :type lat: string
    :return: true if convertion was succesful
    :rtype: bool
    :raises IOError: If the input file doesn't exist.
    """
    try:
        parts = __file__.split(os.sep)
        parts[len(parts) - 1] = "ini_to_json_template.txt"
        temp_adress = os.sep.join(parts)
        template = open(temp_adress,'r')
    except IOError:
        print 'Failed to open template for conversion'
        return False
    else:
        try:
            input_file = open(adr,'r')
        except IOError:
            print 'Failed to open {}'.format(adr)
            return False
        else:
            d_in = json.load(input_file)
            input_file.close()   
            
            output_file = open(res,"w")#clear file
            d_out = read_template(template, d_in, warn, comm, dbg)
            dec = False
            try: 
                dec = d_in['Reaction_module']['Compute_decay']
            except KeyError:
                dec = False
            if(dec == True):
                d_out = decays(d_out, d_in)
            else:
                print "Decay will not be computed." 
            if(mat != ''):
                d_out = mtr.start(d_out, mat, input_file.name)
            
            
            
            new_out = objects_to_arrays( d_out )                        # convert numerical keys to array items
            data = json.dumps(new_out, sort_keys=True, indent=2)        
            data = human_json(data)                                            # humanized json

            #data = array_check(data, alt) #searching for arrays of pairs, also deletion of quation marks
            output_file.write(data)
            #input_file.close()   
            output_file.close()  
    print 'Conversion was sucessful.'
    return True



def read_template(template, d_in, warn, comm, dbg):
    """
    Data are stored in output variable according to
    the template file.

    :param lat: template
    :type lat: array
    :param lat: d_in     input file converted to JSON tree
    :type lat: array
    :param lat: warn
    :type lat: bool
    :param lat: comm
    :type lat: bool
    :param lat: debug
    :type lat: bool
    :return: output variable
    :rtype: array
    """
    d_out = {}
    section = ''
    spaces = ''
    if (warn == True):
        spaces = '  '
    
    transport_on = False
    substances = 2
    
    line_no=0
    for line in template:
        line_no= line_no + 1
        #print 'template line: ', line_no
        line = line.strip()
        if(len(line)<=3):   #empty lines
            continue
        if line.startswith('//'): #comments
            continue
        if line.startswith('['):    #sections
            section = line[1:(line.find(']'))] 
            if(section == 'Settings'):
                break
            continue
        pair = line.split('=')
        pair[0] = pair[0].strip()
        if(len(pair)>1):
            pair[1] = pair[1].strip()
            if pair[1].startswith('NULL'): #useless values
                if (pair[0] == 'N_substances'):#substances exception
                    substances = d_in[section][pair[0]]
                if (dbg == True):
                    #print section + '/' + pair[0] + ' was deleted.'
                    print 'DELETED: ' + spaces + section + '/' + pair[0]
                continue
            if(len(pair)>2): #normal values
                pair[2] = pair[2].strip()
                try: 
                    if (section == 'New keys'):
                        value = ''
                    else:
                        value = d_in[section][pair[0]]          # the value in original file
                except KeyError:
                    if (warn == True):
                        #print section +"/"+pair[0]+" not found."
                        print 'NOT FOUND: ' + section + '/' + pair[0]
                else:
                    if (pair[0] == 'Transport_on'):
                        if (value == True):#transport exception
                            transport_on = True
                            continue
                        else :
                            value = None
                    if (pair[0] == 'Substances'):#substances exception
                        value = value[0:substances]
                    new_adress = pair[1][pair[1].find('"')+1:len(pair[1])]
                    val_type = pair[2][0:pair[2].find('"')]
                    sections = new_adress.split('/')
                    
                    if (section == 'Transport' and transport_on==False
                                and pair[0] != 'Transport_on'):
                        sections[2] += "_save" 
                    if (section == 'New keys' and transport_on==False
                                and pair[0] == 'mobile_p0'):
                        sections[2] += "_save" 
                    # print "error place: ", value, val_type    
                    value = value_type(value, val_type)
                    d_out = update_data(value, sections, d_out)
                    if (dbg == True):
                        #print section + '/' + pair[0] + ' was saved.'
                        print 'SAVED:   ' + spaces + section + '/' + pair[0]
                        print '         ' + spaces + 'new adress: ' + new_adress
                        print '         ' + spaces + 'value: ' + str(value)
                    if(comm == True):
                        try: 
                            value = d_in[section][pair[0] + '_comment']
                        except KeyError:
                            value = ''
                        if(value != ''):
                            sections[len(sections)-1] = 'COMMENT_' + pair[0]
                            d_out = update_data(value, sections, d_out)
    return d_out

    
def objects_to_arrays(d_out):
    """
    convert record with numerical keys to arrays
    recursive
    returns the new tree

    :param lat: d_out     output JSON tree
    :type lat: array
    """
    #print d_out
    
    if (type(d_out) == list):
        new_tree = []
        for item in d_out:
              new_tree.append( objects_to_arrays( item ) )
    
        #print new_tree
        return new_tree
        
    if (type(d_out) == dict):
        
        only_ints=True
        for key, value in d_out.iteritems():
            
            try:
              ret = int( key )
            except ValueError:
              only_ints=False
              break
        
        if ( only_ints ):
            # convert to array
            new_tree=[]
            for key, value in d_out.iteritems():
                index = int( key ) - 1
                assert index >= 0 
                
                while ( len(new_tree) < index + 1  ):
                    new_tree.append(None)
                    
                new_tree[index] = objects_to_arrays( value )
            return new_tree
        else:
            # just recursive calls
            new_tree=d_out
            for key, value in new_tree.iteritems():
                new_tree[key] = objects_to_arrays( value )
            
            return new_tree
            
    else:
        #print new_tree
        return d_out

        

def human_json(json_output_string):
    """
      replace "key" : 
      by key =  
    """
    nicer_keys = re.sub(r'"([a-zA-Z0-9_]*)" *: *',r'\1 = ',json_output_string)
    return nicer_keys
        
        
def decays (res, data):
    """
    Inserts decay data into json file.

    :param lat: res
    :type lat: json
    :param lat: data
    :type lat: json
    :return: updated data
    :rtype: json
    """
    print 'Decays will be computed.'
    try:
        decays_count = data['Reaction_module']['Nr_of_decay_chains']
    except KeyError:
        decays_count = 0
    try:
        for_count = data['Reaction_module']['Nr_of_FoR']
    except KeyError:
        for_count = 0
    
    names_subs = 'Substances'
    names_hl = 'Half_live'
    names_bif = 'Probability'
    names_adress = '/problem/reactions/decays'
    names_kin = 'Kinetic'
    try:
        parts = __file__.split('\\')
        parts[len(parts) - 1] = "ini_to_json_template.txt"
        temp_adress = '\\'.join(parts)
        template = open(temp_adress,'r')
    except IOError:
        print 'Failed to open template for conversion'
        return res
    else :
        setting = 0
        for line in template:
            line = line.strip()
            if (len(line)<=3):   #empty lines
                continue
            if line.startswith('//'): #comments
                continue
            if (setting >= 1):
                value = line[line.find('"')+1:]
                value = value[0:value.find('"')]
                if (setting == 1):
                    names_adress = value
                if (setting == 2):
                    names_subs = value
                if (setting == 3):
                    names_hl = value
                if (setting == 4):
                    names_bif = value
                if (setting == 5):
                    names_kin = value
                setting += 1
            if (line.startswith('[Settings]')):
                setting = 1
    
    current = 1
    for i in range(1, (decays_count + for_count + 1)):
        if (i <= decays_count):
            decay_name = 'Decay_' + str(i)
            iso_count = data[decay_name]['Nr_of_isotopes']
            half_lives = data[decay_name]['Half_lives']
        else :
            decay_name = 'FoReact_' + str(i-decays_count)
            iso_count = 2
            names_hl = names_kin
            half_lives = data[decay_name]['Kinetic_constant']
        if (iso_count == None):
            continue
        try:
            data[decay_name]['Bifurcation_on']
            bifurcation = data[decay_name]['Bifurcation_on']
        except KeyError:
            bifurcation = False
            
        substances = data[decay_name]['Substance_ids']
        
        if(bifurcation == True):
            adress = names_adress + '/' + str(current) + '/' + names_subs
            sections = adress.split('/')
            res = update_data(substances, sections, res)
            
            if (iso_count == 2):
                half = half_lives
            else :
                half = half_lives[0]
            adress = names_adress + '/' + str(current) + '/' + names_hl
            sections = adress.split('/')
            res = update_data(half, sections, res)
            
            bifs = data[decay_name]['Bifurcation']
            adress = names_adress + '/' + str(current) + '/' + names_bif
            sections = adress.split('/')
            res = update_data(bifs, sections, res)
            current += 1
        else :
            for j in range(0, iso_count -1):
                ids = [substances[j], substances[j+1]]
                if (iso_count == 2):
                    half = half_lives
                else :
                    half = half_lives[j]
                adress = names_adress + '/' + str(current) + '/' + names_subs
                sections = adress.split('/')
                res = update_data(ids, sections, res)
                adress = names_adress + '/' + str(current) + '/' + names_hl
                sections = adress.split('/')
                res = update_data(half, sections, res)
                current += 1
    return res
    
    
    
def array_check (data, alt):
    """
    Searches created json for arrays of pairs.
    Also has option to delete quation marks for
    some keys and replace colons with '='

    :param lat: data
    :type lat: json
    :return: updated data
    :rtype: json
    """
    array = False
    array_level = 0
    lines = data.split('\n')
    level = 0
    for i in range(0, len(lines)):
        key = ''
        level = lines[i].find('"')
        if(level==-1):
            level = lines[i].find('}')
        if(level==-1):
            level = lines[i].count(' ')
        line = (lines[i].split('"'))
        if(len(line)>2):
            key = line[1]
            try:
                int(key)
                lines[i] = ' ' * level + '{'
                if(array == False):
                    #print 'opened array on '+str(i)
                    array = True
                    array_level = level
                    lines[i-1] = lines[i-1].replace('{', '[')    
            except ValueError:
                None
            if(alt == True):
                spaces = key.count(' ')
                if(spaces == 0):
                    lines[i] = lines[i].replace('"', "", 2)
                lines[i] = lines[i].replace(':', ' =', 1)
        else:
            if(array == True):
                if(level < array_level):
                    array = False
                    #print 'closed array on '+str(i+1)
                    lines[i] = lines[i].replace('}', ']')   
    data = '\n'.join(lines)
    return data



def value_type (val, val_type):
    """
    Converts value according to specified type.

    :param lat: val
    :type lat: unknown
    :param lat: val_type
    :type lat: string
    :return: value in correct type
    :rtype: unknown
    """
    if(val_type=='double'):
        val = float(val)
    if(val_type=='string'):
        val = str(val)
    if(val_type=='bool'):
        val = bool(val)
    if(val_type=='int'):
        val = int( float(val))
    if(val_type.startswith('enum')):
        parameters = val_type[val_type.find('(')+1:val_type.find(')')]
        pars = parameters.split(',')
        i = 0
        while(i<len(pars)):
            par = pars[i].split('>')
            print "sub:",par[0], par[1], val
            if(str(par[0].strip())==str(val)):
                val = par[1].strip()
                break
            i += 1
    if(val_type.startswith('array of')):
        print "error place: ", val, val_type    
        sub_type=val_type[8:].strip()
        for sub_val in val:
            print "error place: ", sub_val, sub_type        
            sub_val = value_type(sub_val, sub_type)
        
    if(val_type.startswith('data')):
        val = val_type[val_type.find('(')+1:val_type.find(')')]
    return val



def update_data (value, sections, d_out):
    """
    Inserts new data into main output object.

    :param lat: value
    :type lat: unknown
    :param lat: sections
    :type lat: Array
    :param: d_out
    :type: Array
    :return: output data
    :rtype: Array
    """
    len_a = len(sections)
    if(len_a == 2):
        d_out[sections[1]] = value
    if(len_a >= 3):
        if ((sections[1] in d_out) == False):
            d_out[sections[1]] = {} 
        if(len_a==3):
            d_out[sections[1]][sections[2]] = value
    if(len_a >= 4):
        if ((sections[2] in d_out[sections[1]]) == False):
            d_out[sections[1]][sections[2]] = {} 
        if(len_a==4):
            d_out[sections[1]][sections[2]][sections[3]] = value
    if(len_a >= 5):
        if ((sections[3] in d_out[sections[1]][sections[2]]) == False):
            d_out[sections[1]][sections[2]][sections[3]] = {} 
        if(len_a==5):
            d_out[sections[1]][sections[2]][sections[3]][sections[4]] = value
    if(len_a >= 6):
        if ((sections[4] in d_out[sections[1]][sections[2]][sections[3]]) == False):
            d_out[sections[1]][sections[2]][sections[3]][sections[4]] = {} 
        if(len_a==6):
            d_out[sections[1]][sections[2]][sections[3]][sections[4]][sections[5]] = value
    return d_out