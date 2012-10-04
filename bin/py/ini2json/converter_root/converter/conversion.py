# -*- coding: utf-8 -*-
'''
Created on 18.2.2012

@author: Tomáš Košek
'''
level = 1
comment = ""
last_element = 0

def convert (adr, res, comments):
    """
    Starts the convertion of specified file.
    Converted data prints in specified result file.

    :param lat: adr
    :type lat: string
    :param lat: res
    :type lat: string
    :param lat: comments
    :type lat: bool
    :return: true if convertion was succesful
    :rtype: bool
    :raises IOError: If the input file doesn't exist.
    """
    print 'Starting conversion to JSON format.'
    try:
        input_file = open(adr,'r')
    except IOError:
        print 'Failed to open {}'.format(adr)
        return False
    else:
        if(len(res)>5):
            output_file = open(res,"w")
            result_line = ""
            first_line = True
            global level               #number of tabs
            global comment
            global last_element
            output_file.write("{\n")  #starting line of file
            for line in input_file:
                line = line.strip()
                if line.startswith('//'):    #commented
                    result_line = "";
                else:
                    if line.startswith('['):    #header
                        result_line = read_header(line)
                        first_line = True
                    else:   #data
                        if(len(line) > 0):
                            result_line = read_data(line)
                            if(len(result_line) > 0):
                                result_line = tabs() + result_line
                                if (first_line == True):
                                    first_line = False
                                else:
                                    if(comments == False): comment = "";
                                    result_line = ",\n" + result_line + comment
                                    comment = ""
                    output_file.write(result_line)
                    result_line = ""
            if(last_element >= 1):
                level = 1
                result_line = tabs()+"}\n"
                if (first_line == False):
                    result_line = "\n" + result_line
                output_file.write(result_line)
            output_file.write("}")  #end section
        else :
            print "Wrong output file."
            return False
    print 'Conversion to JSON format was sucessful.'
    return True


def read_header (text):
    """
    Reads current line and assembles line for writing.

    :param lat: text
    :type lat: string
    :return: result line in JSON format
    :rtype: string
    """
    result = ""
    global level       
    global comment
    global last_element
    if (last_element >= 1):
        level -= 1
        result = comment + "\n" + tabs() + "},\n\n"
        comment = ""
    result += tabs() + "\""
    result += text[1:text.find(']')] + "\" : {\n"
    level += 1
    last_element = 1
    return result


def read_data (text):
    """
    Reads current line and assembles line for writing.

    :param lat: text
    :type lat: string
    :return: result line in JSON format
    :rtype: string
    """
    pair = text.split("=")
    result = ""
    comment_index = -1
    global comment
    if(len(pair) > 1):
        comment_index = pair[1].find("//")
        #if(len(pair[1][comment_index+2:])==0):    # seems to try catch case of empty comment
        #    comment_index = -1
        #    pair[1] = pair[1][0:comment_index-1]
        if(comment_index>=0):
            #comment = "\t# " + pair[1][comment_index+2:]
            comment = ",\n" + tabs() + '"' + pair[0].strip()+'_comment" : ' + '"' + pair[1][comment_index+2:] + '"'
            pair[1] = pair[1][0:comment_index]
        pair[1] = pair[1].strip()
        # pair[1] = pair[1].replace('""','')
        result = "\"" + pair[0].strip() + "\" : "
        pair[1] = change_words(pair[1])
        pair[1] = pair[1].replace("\t", " ")
        data_type = which_type(pair[1])
        if(data_type == 0):
            if ( len(pair[1])>0 and pair[1][0] == '\"') :
                str_val = pair[1][1:]
                str_val = str_val[0:str_val.find('\"')-1];
            else :
                str_val = pair[1]
            result += "\"" + str_val + "\"" #string
        if(data_type == 1):
            #result += str(float(pair[1])+0) #int + deleting insignificant numbers
            try:
                int(pair[1])
                result += str(int(pair[1])+0)
            except ValueError:
                result += str(float(pair[1])+0)
                
        if((data_type >= 2) and (data_type <= 3)):
            result += pair[1].lower()   #bool and null
        if(data_type == 4):
            result += read_array(pair[1])
    return result


def read_array (text):
    """
    Reads current line and assembles line for writing.

    :param lat: text
    :type lat: string
    :return: result line in JSON format
    :rtype: string
    """
    result = "["
    field = text.split(" ")
    first = True
    for value in field:
        if (len(value) > 0):
            if (first == True):
                first = False
            else :
                result += ", "
            data_type = which_type(value)
            if(data_type == 0):
                result += "\"" + value + "\"" #string
            if(data_type == 1):
                #result += str(float(value)+0) #int + deleting insignificant numbers
                try:
                    int(value)
                    result += str(int(value)+0)
                except ValueError:
                    result += str(float(value)+0)
            if((data_type >= 2) and (data_type <= 3)):
                result += value.lower()   #bool and null
    result += "]"
    return result


def change_words(word):
    """
    Checks for several specified words in value field.
    On, yes and off, no are changed to true and false. 

    :param lat: word
    :type lat: string
    :return: corrected word
    :rtype: string
    """
    word = word.strip()
    if (word.isalnum()):
        if (word.lower() == "on"):
            word = "true"
        if (word.lower() == "yes"):
            word = "true"
        if (word.lower() == "off"):
            word = "false"
        if (word.lower() == "no"):
            word = "false"
    return word


def which_type (value):
    """
    Reads current value and searches which data format is used.

    :param lat: value
    :type lat: string
    :return: number of current type
    :rtype: int
    """
    try:
        float(value)
        return 1 #float or int
    except ValueError:
        if ( len(value)>0 and value[0] == '\"') :
            return 0
        if((value == "true") or (value == "false")):
            return 2 #boolean
        if(value.lower() == "null"):
            return 3 #null
        spaces = value.split(" ")
        if (len(spaces) > 1):
            return 4
        return 0


def tabs ():
    """
    Simply adds 'level' tabs for better readability

    :return: string with tabs characters
    :rtype: string
    """
    global level
    return (" ")