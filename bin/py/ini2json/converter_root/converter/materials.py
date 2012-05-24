# -*- coding: utf-8 -*-
'''
Created on 15.5.2012

@author: Tomáš Košek
'''

import completion as comp

def start (data, adress, input_adress):
    """
    Starts the convertion of material file.
    Result data are sent back for final operations.

    :param lat: data
    :type lat: array
    :param lat: adress
    :type lat: string
    :param lat: input_adress
    :type lat: string
    :return: data with addition of values from materials
    :rtype: array
    :raises IOError: If the input file doesn't exist.
    """
    if(adress == '*'):
        try: 
            adress = data['material']
        except KeyError:
            print 'Material file not found.'
            return data
        temp = input_adress.split('/')
        if (adress.startswith('./')):
            adress = adress[2:]
        temp[len(temp)-1] = adress
        adress = '/'.join(temp)
    
    try:
        input_file = open(adress,'r')
    except IOError:
        print 'Failed to open {}'.format(adress)
        return data
    else:
        section = ''
        mats = {}
        stor = {}
        geom = {}
        sorp = {}
        frac = {}
        dual = {}
        
        for line in input_file:
            if(line.startswith('#')):
                continue
            line = line.replace('\t',' ')
            line = line.strip()
            while (line.count('  ')>0):
                line = line.replace('  ',' ')
            
            if(line.startswith('$End')):
                section = ''
                continue
            elif (line.startswith('$')):
                section = line[1:]
                counter = 1
                continue
            if((section != '') and (len(line)>0)):
                line_data = line.strip().split(' ')
                for i in range(0, len(line_data)):
                    try:
                        line_data[i] = int(line_data[i].strip())
                    except ValueError:
                        line_data[i] = float(line_data[i].strip())
                
                if(section == 'Materials'):
                    if(len(line_data) > 2):
                        mats[counter] = sec_materials(line_data)
                    else:
                        continue
                    
                if(section == 'Storativity'):
                    stor[counter] = {'material' : line_data[0], 
                                     'value' : line_data[1]}
                    
                if(section == 'Geometry'):
                    geom[counter] = {'material' : line_data[0], 
                                     'value' : line_data[2]}
                    
                if(section == 'Sorption'):
                    sorp[counter] = sec_sorption(line_data)
                    
                if(section == 'SorptionFraction'):
                    frac[counter] = {'material' : line_data[0], 
                                     'value' : line_data[1]}
                
                if(section == 'DualPorosity'):
                    dual[counter] = sec_dual(line_data)
                    
                counter += 1
    
    data['problem']['primary_equation']['conductivity'] = mats
    data['problem']['primary_equation']['storativity'] = stor
    data['problem']['primary_equation']['cross_area'] = geom
    data['problem']['secondary_equation']['sorption'] = sorp
    data['problem']['secondary_equation']['sorption_fraction'] = frac
    data['problem']['secondary_equation']['dual_porosity'] = dual
    return data
            
def sec_materials (values):
    """
    This returns correct form of data in materials section.

    :param lat: values
    :type lat: array
    :return: array with values from section materials
    :rtype: array
    """
    ret = {}
    ret['material'] = values[0]
    if(values[1] == 11):
        ret['value'] = values[2]
    if(values[1] == 21):
        ret['value'] = [values[2], values[2], 0]
    if(values[1] == 22):
        ret['value'] = [values[2], values[3], 0]
    if(values[1] == 23):
        ret['value'] = [values[2], values[3], values[4]]
    if(values[1] == 31):
        ret['value'] = [[values[2], 0, 0], 
                        [0, values[2], 0], 
                        [0, 0, values[2]]]   
    if(values[1] == 33):
        ret['value'] = [[values[2], 0, 0], 
                        [0, values[3], 0], 
                        [0, 0, values[4]]]
    if(values[1] == 36):
        ret['value'] = [[values[2], values[5], values[6]],
                        [values[5], values[3], values[7]],
                        [values[6], values[7], values[4]]]
    return ret
            
            
def sec_dual (values):
    """
    This returns correct form of data in dual porosity section.

    :param lat: values
    :type lat: array
    :return: array with values from section materials
    :rtype: array
    """
    ret = {}
    ret['material'] = values[0]
    ret['mobile'] = values[1]
    ret['immobile'] = values[2]
    ret['coef'] = [values[i] for i in range(3, len(values))]
    return ret   


def sec_sorption (values):
    """
    This returns correct form of data in sorption section.

    :param lat: values
    :type lat: array
    :return: array with values from section materials
    :rtype: array
    """
    ret = {}
    ret['material'] = values[0]
    ret['substance'] = values[1]
    enum = 'enum(1>Equilibrium, 2>Freundlich, 3>Langmuir)'
    ret['type'] = comp.value_type(values[2], enum)
    coefs = 1
    if(values[2]>1):
        coefs = 2
    ret['coef'] = [values[i] for i in range(3, 3+coefs)]
    return ret         