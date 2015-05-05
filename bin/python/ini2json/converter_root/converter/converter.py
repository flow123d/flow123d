# -*- coding: utf-8 -*-
'''
Created on 8.2.2012

@author: Tomáš Košek
'''
import conversion as cnv
import completion as comp

import argparse

def main():

    #parsing
    print 'Program started.'
    parser = argparse.ArgumentParser(description=
                                     'Enter parameters for conversion.')
    parser.add_argument('-ini', action='store', default='', dest='file_ini',
                        help='Adress or name of source file in ini.')
    parser.add_argument('-json', action='store', default='', dest='file_json',
                        help='Adress or name of source file in json. (UNUSED)')
    parser.add_argument('-final', action='store', default='', dest='file_final',
                        help='Adress or name of converted file.')
    parser.add_argument('-mat', action='store', default='', dest='file_material',
                        help='Adress or name of material file. Type * for automatic file detection.')
    parser.add_argument('-sw', action='store_true', dest='warnings',
                        help='Show not found messages during conversion.')
    parser.add_argument('-nc', action='store_false', dest='comments',
                        help='Disable saving of comments.')
    parser.add_argument('-db', action='store_true', dest='debug',
                        help='Show saved values (used for debugging).')
    parser.add_argument('-nd', action='store_false', dest='altering',
                        help='Disable alterations to json file (quotation marks, colons).')
    
    args = vars(parser.parse_args())
    #print args
    
    #checking for empty parameters
    f1 = args['file_ini']
    f2 = args['file_json']
    f3 = args['file_final']
    converted = False
    if (f1 == ''):
        if (f2 == ''):
            print 'Error: no file name inserted.'
        else :
            converted = True
            if (f3 == ''):
                f3 = f2[0:-3] + 'json'
                print 'Name of final file: ' + f3
    else :
        if (f2 == ''):
            f2 = f1[0:-4] + '_temp.json'
            print 'Name of source file in json: ' + f2
        if (f3 == ''):
            f3 = f1[0:-3] + 'json'
            print 'Name of final file: ' + f3
        converted = cnv.convert(f1, f2, args['comments'])   #converts inifile into JSON
    
    #alter json file
    if (converted == True):
        print 'Starting final conversion.'
        complete = comp.complete(f2, f3, args['warnings'], args['comments'], 
                                 args['debug'], args['altering'],
                                 args['file_material'])
    print 'Program terminated.'
    pass

if __name__ == '__main__':
    main()