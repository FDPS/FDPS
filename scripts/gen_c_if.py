#!/usr/bin/env python
# -*- coding: utf-8 -*-
#=================================
#   Module import
#=================================
try:
    import os
    import sys
    import re
    import struct
    import math
except ImportError:
    print("Modules os,sys,struct,math,string are not found.")
    quit()

try:
    import string
    import random
except ImportError:
    print("Modules string, random are not found.")
    quit()

try:
    import shutil
except ImportError:
    print("Module shutil is not found.")
    quit()


try:
    import argparse 
except ImportError:
    print("Module argparse is not found.")
    quit()

try:
    import itertools
except ImportError:
    print("Module itertools is not found.")
    quit()

try:
    import collections
except ImportError:
    print("Module collections is not found.")
    quit()

try:
    import textwrap
except ImportError:
    print("Module textwrap is not found.")
    quit()

try:
    import subprocess
except ImportError:
    print("Module subprocess is not found.")
    quit()

#================================
#   Class definition
#================================
class File_Data:
    def __init__(self):
        self.name    = ""
        self.content = []
        self.structures = []

class Structure_Data:
    def __init__(self):
        # line numbers containing various keywords used when defining
        # a C structure 
        self.line_num_struct_key = -1 # the line number containing `struct`
        self.line_num_right_curly_brace = -1 # the line number containing `}`
        self.line_num_left_curly_brace = -1 # the line number containing `{`
        self.line_num_tag_name = -1 # the line number containing tag name

        # The start and end points of `struct` in the relevant line
        self.start_struct_key = -1
        self.end_struct_key = -1

        # The start points of `{` and `}` in the relevant lines
        self.start_left_curly_brace = -1
        self.start_right_curly_brace = -1

        # The start and end points of tag name in the relevant line
        self.start_tag_name = -1
        self.end_tag_name = -1
        self.tag_name = ""

        # The content of FDPS directive that specifies that this C-structure
        # corresponds to which user-defined type.
        self.fdps_str_dir = ""

        # The list of FDPS directives that specifies how data copies are
        # performed in FDPS.
        self.fdps_meth_dirs = []

        # Variable declaration parts
        self.var_decl_parts = []

        # Member variables and FDPS directive for member variables
        self.members      = []

class Variable_Declaration_Data:
    def __init__(self):
        self.line_num_start = -1 # line number containing the 1st [a-zA-Z_] character
        self.pos_start      = -1 # position of the 1st [a-zA-Z_] character in the relevant line
        self.line_num_end   = -1 # line number contaning the character `;`.
        self.pos_end        = -1 # position of the character `;` in the relavant line
        self.text           = "" # text data of the variable desclaration part

class Member_Data:
    def __init__(self):
        self.name       = ""
        self.attrib     = ""
        self.line_num   = -1
        self.data_type  = ""
        self.is_array   = False
        self.array_dim  = -1

class Automatic_Generator:
    def __init__(self,input_files,output_dir,dim_num):
        #=======================================================================
        # Display copyright
        text = """
        ################################################

        Python script to generate FDPS C interface 
        (C) Copyright 2018  FDPS developer team. 

        Log information 
           Python version   = {VERSION}
           Input files      = {INPUF_FILES}
           Output directory = {OUTPUT_DIR}
           Dimension number = {DIM_NUM}

        ################################################
        """.format(VERSION=sys.version_info, \
                   INPUF_FILES=input_files, \
                   OUTPUT_DIR=output_dir, \
                   DIM_NUM=dim_num)
        text = text[1:].rstrip() + "\n"
        text = textwrap.dedent(text)
        print(text)

        #=======================================================================
        # Set the dimension number of the space
        self.__DIM_NUM = dim_num

        #=======================================================================
        # Data-type conversion dictionary, and 
        # The list of data types usable in Fortran codes
        # (For details, see https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING)
        # 
        # * List of interoperable types
        #-----------------------------------------------------------------------
        # integer(kind=C_INT)                   int 
        # integer(kind=C_SHORT)                 short int 
        # integer(kind=C_LONG)                  long int 
        # integer(kind=C_LONG_LONG)             long long int 
        # integer(kind=C_SIGNED_CHAR)           signed char/unsigned char 
        # integer(kind=C_SIZE_T)                size_t 
        # integer(kind=C_INT8_T)                int8_t 
        # integer(kind=C_INT16_T)               int16_t 
        # integer(kind=C_INT32_T)               int32_t 
        # integer(kind=C_INT64_T)               int64_t 
        # integer(kind=C_INT_LEAST8_T)          int_least8_t 
        # integer(kind=C_INT_LEAST16_T)         int_least16_t 
        # integer(kind=C_INT_LEAST32_T)         int_least32_t 
        # integer(kind=C_INT_LEAST64_T)         int_least64_t 
        # integer(kind=C_INT_FAST8_T)           int_fast8_t 
        # integer(kind=C_INT_FAST16_T)          int_fast16_t 
        # integer(kind=C_INT_FAST32_T)          int_fast32_t 
        # integer(kind=C_INT_FAST64_T)          int_fast64_t 
        # integer(kind=C_INTMAX_T)              intmax_t 
        # integer(kind=C_INTPTR_T)              intptr_t 
        # real(kind=C_FLOAT)                    float 
        # real(kind=C_DOUBLE)                   double 
        # real(kind=C_LONG_DOUBLE)              long double 
        # complex(kind=C_FLOAT_COMPLEX)         float _Complex 
        # complex(kind=C_DOUBLE_COMPLEX)        double _Complex 
        # complex(kind=C_LONG_DOUBLE_COMPLEX)   long double _Complex 
        # logical(kind=C_BOOL)                  _Bool 
        # character(kind=C_CHAR)                char 
        #-----------------------------------------------------------------------
        self.__cDT_to_ftnDT = {"int":"integer(kind=c_int)", \
                               "short":"integer(kind=c_short)", \
                               "short int":"integer(kind=c_short)", \
                               "int short":"integer(kind=c_short)", \
                               "long":"integer(kind=c_long)", \
                               "long int":"integer(kind=c_long)", \
                               "int long":"integer(kind=c_long)", \
                               "long long":"integer(kind=c_long_long)", \
                               "long long int":"integer(kind=c_long_long)", \
                               "fdps_s32":"integer(kind=c_int)", \
                               "fdps_s64":"integer(kind=c_long_long)", \
                               "signed char":"integer(kind=c_signed_char)", \
                               "size_t":"integer(kind=c_size_t)", \
                               "int8_t":"integer(kind=c_int8_t)", \
                               "int16_t":"integer(kind=c_int16_t)", \
                               "int32_t":"integer(kind=c_int32_t)", \
                               "int64_t":"integer(kind=c_int64_t)", \
                               "int_least8_t":"integer(kind=c_int_least8_t)", \
                               "int_least16_t":"integer(kind=c_int_least16_t)", \
                               "int_least32_t":"integer(kind=c_int_least32_t)", \
                               "int_least64_t":"integer(kind=c_int_least64_t)", \
                               "int_fast8_t":"integer(kind=c_int_fast8_t)", \
                               "int_fast16_t":"integer(kind=c_int_fast16_t)", \
                               "int_fast32_t":"integer(kind=c_int_fast32_t)", \
                               "int_fast64_t":"integer(kind=c_int_fast64_t)", \
                               "intmax_t":"integer(kind=c_intmax_t)", \
                               "intptr_t":"integer(kind=c_intptr_t)", \
                               "float":"real(kind=c_float)", \
                               "double":"real(kind=c_double)", \
                               "long double":"real(kind=c_long_double)", \
                               "fdps_f32":"real(kind=c_float)", \
                               "fdps_f64":"real(kind=c_double)", \
                               "float _Complex":"complex(kind=c_float_complex)", \
                               "double _Complex":"complex(kind=c_double_complex)", \
                               "long double _Complex":"complex(kind=c_long_double_complex)", \
                               "_Bool":"logical(kind=c_bool)", \
                               "bool":"logical(kind=c_bool)", \
                               "char":"character(kind=c_char)", \
                               "fdps_f32vec":"type(fdps_f32vec)", \
                               "fdps_f64vec":"type(fdps_f64vec)", \
                               "fdps_f32mat":"type(fdps_f32mat)", \
                               "fdps_f64mat":"type(fdps_f64mat)"}
        self.__integer_cDT = frozenset(["int", \
                                        "short", \
                                        "short int", \
                                        "int short", \
                                        "long", \
                                        "long int", \
                                        "int long", \
                                        "long long", \
                                        "long long int", \
                                        "fdps_s32", \
                                        "fdps_s64", \
                                        "signed char", \
                                        "size_t", \
                                        "int8_t", \
                                        "int16_t", \
                                        "int32_t", \
                                        "int64_t", \
                                        "int_least8_t", \
                                        "int_least16_t", \
                                        "int_least32_t", \
                                        "int_least64_t", \
                                        "int_fast8_t", \
                                        "int_fast16_t", \
                                        "int_fast32_t", \
                                        "int_fast64_t", \
                                        "intmax_t", \
                                        "intptr_t"])
        self.__float_cDT   = frozenset(["float", \
                                        "double", \
                                        "long double", \
                                        "fdps_s32", \
                                        "fdps_f64", \
                                        "float _Complex", \
                                        "double _Complex", \
                                        "long double _Complex", \
                                        "fdps_f32vec", \
                                        "fdps_f64vec", \
                                        "fdps_f32mat", \
                                        "fdps_f64mat"])
        self.__boolean_cDT = frozenset(["_Bool", \
                                        "bool"])
        self.__character_cDT = frozenset(["char"])
        self.__usable_data_types = frozenset(self.__cDT_to_ftnDT.keys())
        self.__reserved_words_in_c = set()
        for data_type in self.__usable_data_types:
            items = data_type.split()
            for item in items:
                self.__reserved_words_in_c.add(item)
        #=======================================================================
        # Usable FDPS directive keywords
        self.__usable_fdps_str_dirs = frozenset(["FP", \
                                                 "EPI", \
                                                 "EPJ", \
                                                 "Force"])
        self.__usable_fdps_mbr_dirs = frozenset(["id", \
                                                 "position", \
                                                 "velocity", \
                                                 "charge", \
                                                 "rsearch"])
        self.__usable_fdps_meth_dirs = frozenset(["copyFromForce", \
                                                  "copyFromForcePM", \
                                                  "copyFromFP", \
                                                  "clear"])
        #=======================================================================
        # C file manager
        self.files      = []

        #=======================================================================
        # List of FDPS user defined types
        self.__FPs = []
        self.__EPIs = []
        self.__EPJs = []
        self.__Forces = []

    def __print_warning(self,msg):
        output_msg = "\033[33;1m" + "[warning] " + msg + "\033[m" # yellow
        print(output_msg)

    def __print_error(self,msg):
        output_msg = "\033[31;1m" + "[error] " + msg + "\033[m" # red
        print(output_msg)

    def __print_checkpoint(self,msg):
        output_msg = "\033[32;1m" + "[check point] " + msg + "\033[m" # green
        print(output_msg)

    def __print_dbgmsg(self,msg):
        output_msg = "\033[35;1m" + "[debug] " + msg + "\033[m" # magenta
        print(output_msg)

    def __get_random_string(self,n):
        c = string.ascii_lowercase + string.ascii_uppercase + string.digits
        return ''.join([random.choice(c) for i in range(n)])

    def __get_header_cmt(self,text):
        # This function checks if `*/` exists or not. And if it exists,
        # this function divides a given text into a divided comment part (header)
        # and the remaining part (body). If cannot find `*/`, it means that
        # all of this line is a divided comment.
        closer_is_founded = 0
        header = ""
        body   = ""
        m = re.search(r"^(?P<header>.*?)\*/",text)
        if (m):
            closer_is_founded = 1
            header = m.string[m.start():m.end()]
            body   = m.string[m.end():]
        else:
            header = text
            body   = ""
        # Check if header contains FDPS directives
        m = re.search(r"\$fdps",header)
        if (m):
            return header.strip(), body
        else:
            if (closer_is_founded):
                return "*/", body
            else:
                return "", body

    def __get_pos_left_open_cmt(self,text):
        # This function checks if there is a left open comment in the given text.
        # If it is found, the function returns the position of `*` of '*/'.
        # Otherwise, it returns -1.
        pos = -1
        if (len(text)):
            is_in_cmt = 0
            char_prevs = [" "," "] #
            for i in range(0,len(text)):
                char = text[i]
                str_tgt = char_prevs[1] + char
                if (is_in_cmt == 0):
                    if (str_tgt == "/*"):
                        is_in_cmt = 1
                    elif (str_tgt == "*/"):
                        pos = i - 1
                        break
                    elif (str_tgt == "//"):
                        # In this case, the remaining text must be a comment.
                        # Therefore, we do not need to continue the check.
                        break
                else:
                    if (str_tgt == "*/"):
                        is_in_cmt = 0
                        char_prevs = [" "," "] # reset
                        # This reset is required to prevent a wrong detection of
                        # double slash comment. Here is an example:
                        # (ex1) z = x  /* aaa *// y;
                        # The last slash is not a part of a double slash comment.
                        # In the case above, the last slash represents division symbol.
                        # If we reset `char_prevs` when the control reaches the second
                        # slash, we can prevent this wrong detection.
                        # Confusingly, a double slash in the following code snippet
                        # does represent a C-comment:
                        # (ex2) z = x *// y
                        #           y;
                char_prevs[0] = char_prevs[1]
                char_prevs[1] = char
        return pos 

    def __find_pos_cmt2(self,text):
        # This function searches /**/ comments in a given text
        # and return a list of the index ranges (start,end) of
        # the detected /**/ comment (the 1st character in the 
        # line corresponds to 0).
        # This function assumes that the head of the given text
        # is not a part of other /**/ comment.
        rslt = []
        if (len(text)):
            pos = 0
            while 1:
                m = re.search(r"/\*(?P<cmt>.*?)\*/",text[pos:])
                if (m):
                    start = pos + m.start()
                    end = pos + m.end() - 1
                    rslt.append((start,end))
                    pos = end + 1
                else:
                    break
            # At last, we check if there is a right open comment.
            m = re.search(r"/\*(?P<cmt>.*?)$",text[pos:])
            if (m):
                start = pos + m.start()
                end = pos + m.end() - 1
                rslt.append((start,end))
        return rslt

    def __find_pos_non_cmt(self,text,start_srch=0,end_srch=-1):
        # This function returns a list of position ranges of text
        # where strings of characters are not C comments.
        # Error handling
        len_text = len(text)
        if (end_srch == -1):
            end_srch = len_text - 1
        if ((start_srch < 0) or (end_srch > len_text-1)):
            msg = "The values of the default arguments of __find_pos_non_cmt are invalid!"
            self.__print_error(msg)
            sys.exit()
        # Examine the parts of the text which are not C comment.
        pos_cmt = self.__find_pos_cmt2(text)
        pos_non_cmt = []
        if pos_cmt:
            # First, find parts of the text which are not C comments.
            start = 0
            pos_non_cmt_cand = []
            for start_cmt,end_cmt in pos_cmt:
                end = start_cmt - 1
                pos_non_cmt_cand.append([start,end])
                start = end_cmt + 1
            end = len(text) - 1
            pos_non_cmt_cand.append([start,end])
            # Then, discord offset
            while pos_non_cmt_cand:
                start, end = pos_non_cmt_cand.pop()
                if (end < start_srch):
                    continue
                elif ((start <= start_srch) and (start_srch <= end)):
                    start = start_srch
                if (end_srch < start):
                    continue
                elif ((start <= end_srch) and (end_srch <= end)):
                    end = end_srch
                pos_non_cmt.append([start,end]) 
        else:
            pos_non_cmt.append([start_srch,end_srch])
        return pos_non_cmt    

    def __change_cmt_type(self,text):
        # First, calculate index ranges in which strings of characters
        # may not be C-comment.
        pos_non_cmt = []
        pos_cmt = self.__find_pos_cmt2(text)
        if pos_cmt:
            start = 0
            while pos_cmt:
                start_cmt, end_cmt = pos_cmt.pop()
                end = start_cmt - 1
                if (start <= end):
                    pos_non_cmt.append((start,end))
                start = end_cmt + 1
        else:
            pos_non_cmt.append((0,len(text)-1))
        # Next, search `//` in the probably-non-comment areas of the text
        indx = -1 # the position of the 2nd slash
        if pos_non_cmt:
            while pos_non_cmt:
                flag_break = False
                start, end = pos_non_cmt.pop()
                char_prev = ""
                for i in range(start,end+1):
                    char = text[i]
                    str_tgt = char_prev + char
                    if (str_tgt == "//"):
                        indx = i 
                        flag_break = True
                        break
                    char_prev = char
                if (flag_break):
                    break
        # Finally, we replace // comment to /**/ comment
        if (indx != -1):
            return text[0:indx] + "*" + text[indx+1:] + "*/"
        else:
            return text

    def __delete_unnecessary_cmt(self,text):
        # Here, we check the contents of the closed /**/ comments in this text.
        rslt = ""
        if (len(text)):
            pos = 0
            while 1:
                m = re.search(r"/\*(?P<cmt>.*?)\*/",text[pos:])
                if (m):
                    cmt = m.group('cmt')
                    start = pos + m.start()
                    end = pos + m.end()
                    m = re.search(r"\$fdps",cmt)
                    if (m):
                        rslt += text[pos:end]
                    else:
                        # In this case, we delete this comment
                        rslt += text[pos:start]
                    pos = end
                else:
                    rslt += text[pos:]
                    break
        return rslt

    def __get_footer_cmt(self,text):
        # This function checks if a divided /**/ comment exists at the end of text.
        # If it exists, this function divides text into the leading part and 
        # the divided comment.
        pos_footer_cmt = -1
        is_in_cmt = 0
        pos = 0
        char_prev = ""
        for char in text:
            str_tgt = char_prev + char
            if (is_in_cmt == 0):
                if (str_tgt == "/*"):
                    is_in_cmt = 1
                    pos_footer_cmt = pos - 1
            else:
                if (str_tgt == "*/"):
                    is_in_cmt = 0
                    pos_footer_cmt = -1
            char_prev = char
            pos += 1
        # Set body, footer 
        if (pos_footer_cmt != -1):
            body   = text[:pos_footer_cmt]
            footer = text[pos_footer_cmt:]
        else:
            body   = text
            footer = ""
        return body.strip(), footer.strip()

    def __find_regex(self,text,regex,offset=0,\
                     check_brace_pairing=False, \
                     num_left_brace_ini=0, \
                     num_right_brace_ini=0, \
                     debug_mode=False):
        # This function seeks the 1st string of characters which matches
        # the given regular expression and which are not C-comment.
        # If a string of characters is found, the function returns
        # the position information about the string, i.e., (start,end).
        # If check_brace_pairing is True, the function counts the numbers
        # of braces, `(` and ')', in the part of the text from offset to
        # just before the matched string and returns them. If the matched
        # string is not found, the numbers of braces in the part of
        # text[offset:] are returned.
        
        # Error handling
        len_text_org = len(text)
        if (offset > len_text_org):
            errmsg = "A invalid value is specified for offset in __find_regex.\n" \
                   + "text   = {0}\n".format(text) \
                   + "offset = {0}\n".format(offset)
            self.__print_dbgmsg(errmsg)
            sys.exit()
        # Initialize
        num_left_brace = num_left_brace_ini
        num_right_brace = num_right_brace_ini
        # Main body of the funtion
        start_tgt = -1
        end_tgt   = -1
        if (len_text_org):
            # We first get the ranges of the positions of the non-comment parts.
            pos_non_cmt = self.__find_pos_non_cmt(text,offset)
            # Then, we search each part.
            while pos_non_cmt:
                start_srch, end_srch = pos_non_cmt.pop()
                flag_break = False
                while start_srch <= end_srch:
                    m = re.search(regex,text[start_srch:end_srch+1])
                    if (m):
                        start_tgt = start_srch + m.start('pat')
                        end_tgt   = start_srch + m.end('pat') - 1
                        flag_break = True
                        if (check_brace_pairing):
                            # In this case, we counts # of `(` and ')' in the
                            # text before the initial letter of the matched string.
                            # And if # of '(' and ')' are not the same, we discord
                            # the matched string.
                            for i in range(start_srch,start_tgt):
                                char = text[i]
                                if (char == "("):
                                    num_left_brace += 1
                                elif (char == ")"):
                                    num_right_brace += 1
                            if (num_left_brace != num_right_brace):
                                flag_break = False
                                start_srch = end_tgt + 1
                                start_tgt = -1
                                end_tgt = -1
                            else:
                                break # escape from the inner while-loop
                    else:
                        # In this case, we cannot find a matched string
                        # in this range of the text. So, we move to
                        # the next range.
                        if (check_brace_pairing):
                            for i in range(start_srch,end_srch+1):
                                char = text[i]
                                if (char == "("):
                                    num_left_brace += 1
                                elif (char == ")"):
                                    num_right_brace += 1
                        break # escape from the inner while-loop
                    if (flag_break):
                        break
        return start_tgt,end_tgt,num_left_brace,num_right_brace

    def __find_struct_keyword(self,text,       \
                              num_left_brace,  \
                              num_right_brace, \
                              offset=0):
        # This function seeks the keyword `struct` that is not comment
        # in a given text. If the keyword is found, the function returns
        # the positional information about the keyword.
        regex = r"\s*(?P<pat>struct)\s*" 
        return self.__find_regex(text,regex,offset,True,\
                                 num_left_brace,num_right_brace)

    def __find_tag_name(self,text,offset=0):
        # This function seeks the 1st string of characters which matches
        # the regular expression [a-zA-z_][a-zA-Z0-9_]* and which are not
        # comment. If a string of characters is found, the function returns
        # the positional information about the string.
        regex = r"(?P<pat>[a-zA-Z_][a-zA-Z0-9_]*)"
        start, end, n_lbrace, n_rbrace = self.__find_regex(text,regex,offset,debug_mode=True)
        return start, end

    def __get_fdps_str_dir(self,text):
        # This function first checks if a given text satisfies the format
        # of FDPS directive for structure. If it is the valid directive,
        # the function extracts it and returns it.
        # Note that this function assumes that comment mark /**/ are
        # eliminated in advance.
        fdps_dir = ""
        m = re.match(r"\$fdps\s+(?P<keyword_list>\S+)",text.strip())
        if (m):
            keyword_list = m.group('keyword_list')
            keys = keyword_list.split(',')
            for key in keys:
                key = key.strip()
                if key not in self.__usable_fdps_str_dirs:
                    msg = "An unknown keyword is used in a FDPS directive!\n" \
                        + "The keyword used is {0}".format(key)
                    self.__print_error(errmsg)
                    sys.exit()
            fdps_dir = m.group()
        else:
            msg = "There are no keywords in this FDPS directive!\n" \
                + "Please check the following part of the code:\n" \
                + "{0}".format(text.strip())
            self.__print_error(msg)
        return fdps_dir

    def __is_fdps_str_dir(self,text):
        # This function checks if a given text satisfies the format
        # of FDPS directive for structure. The result is returned
        # as a boolean value.
        # Note that this function assumes that comment mark /**/ are
        # eliminated in advance.
        m = re.match(r"\$fdps\s+(?P<keyword_list>\S+)",text.strip())
        if (m):
            keyword_list = m.group('keyword_list')
            keys = keyword_list.split(',')
            for key in keys:
                key = key.strip()
                if key not in self.__usable_fdps_str_dirs:
                    return False
            return True
        else:
            return False

    def __get_fdps_meth_dir_cand(self,text):
        # This function first checks if a given text satisfies the format
        # of FDPS directive for method IN A RELATIVELY SIMPLE WAY.
        # If the text passes this simple check, the function extracts
        # a FDPS directive candidate and returns it.
        # Note that this function assumes that comment mark /**/ are
        # eliminated in advance.
        m = re.match(r"\$fdps\s+(?P<meth_name>[a-zA-Z]+)\s.*",text.strip())
        if (m):
            meth_name = m.group('meth_name').strip()
            if meth_name in self.__usable_fdps_meth_dirs:
                return m.group()
            else:
                return ""
        else:
            return ""

    def __get_fdps_mbr_dir(self,text):
        # This function checks if a given text satisfies the format
        # of FDPS directive for member variable.
        # If the text passes the check, the function extracts 
        # keyword such as id, position, etc. and return it.
        # Otherwise, return empty string of characters.
        # Note that this function assumes that comment mark /**/ are
        # eliminated in advance.
        m = re.match(r"^\$fdps\s+(?P<key>[a-z]+)$",text.strip())
        if (m):
            key = m.group('key')
            if key in self.__usable_fdps_mbr_dirs:
                return key
            else:
                return ""
                # In this case, this directive probably is a FDPS directive
                # that specifies the way of data processing in FDPS such as
                # $fdps clear. So, we ignore this case.
        else:
            return ""

    def __struct_existence_check(self,tag_name):
        bool_val = False
        for f in self.files:
            for s in f.structures:
                if (tag_name == s.tag_name):
                    bool_val = True
        return bool_val 

    def __member_existence_check(self,tag_name,member_name):
        bool_val = False
        for f in self.files:
            for s in f.structures:
                if (tag_name == s.tag_name):
                    for mbr in s.members:
                        if (member_name == mbr.name):
                            bool_val = True
        return bool_val

    def __is_s32(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "integer(kind=c_int)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_s64(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "integer(kind=c_long_long)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f32(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "real(kind=c_float)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f64(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "real(kind=c_double)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f32vec(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "type(fdps_f32vec)"))
        conditions.append((ftnDT == "type(kind=c_float)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        return any(conditions)

    def __is_f64vec(self,mbr):
        ftnDT = self.__cDT_to_ftnDT[mbr.data_type]
        conditions = []
        conditions.append((ftnDT == "type(fdps_f64vec)"))
        conditions.append((ftnDT == "type(kind=c_double)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        return any(conditions)

    def __is_integer(self,data_type):
        if data_type in self.__integer_cDT:
            return True
        else:
            return False

    def __is_float(self,data_type):
        if data_type in self.__float_cDT:
            return True
        else:
            return False

    def __is_boolean(self,data_type):
        if data_type in self.__boolean_cDT:
            return True
        else:
            return False

    def __is_character(self,data_type):
        if data_type in self.__character_cDT:
            return True
        else:
            return False

    def __get_data_type(self,tag_name,mbr_name):
        for f in self.files:
            for s in f.structures:
                if (tag_name == s.tag_name):
                    for mbr in s.members:
                        if (mbr_name == mbr.name):
                            return mbr.data_type
        return ""

    def __get_fdpsDT(self,mbr):
        # This method converts C data type to Fortran data type used in FDPS.
        # This method works correctly only if data in `mbr` is complete.
        if (mbr.is_array):
            if (mbr.attrib not in self.__usable_fdps_mbr_dirs):
                # normal array cases
                data_type = self.__cDT_to_ftnDT[mbr.data_type]
                array_dim = mbr.array_dim
            else:
                # `position` or `velocity`
                if (self.__is_f32vec(mbr)):
                    data_type = "type(fdps_f32vec)"
                    array_dim = 0
                else:
                    data_type = "type(fdps_f64vec)"
                    array_dim = 0
        else:
            data_type = self.__cDT_to_ftnDT[mbr.data_type]
            array_dim = 0
        return data_type,array_dim

    def __add_indent(self,input_text,indent):
        output_text = ""
        items = input_text.split('\n')
        item_num = 0
        for item in items:
            item_to_be_added = (indent + item + "\n")
            if (item_num < len(items)):
                output_text += item_to_be_added
            else:
                if (item != ''):
                    output_text += item_to_be_added
        return output_text

    def read_files(self,file_names):
        # Output to stdout
        self.__print_checkpoint("reading user's C files")
        for fname in file_names:
            print("--- reading {0}".format(fname))
            # Open file
            fp = open(fname,'r')
            # Create a file object
            self.files.append(File_Data())
            # Get the reference for the current file
            f = self.files[-1]
            # Set the file name
            f.name = fname
            # Read the file
            is_in_cmt = 0
            line_buf  = ""
            for line in fp:
                # Delete LF
                line = line.rstrip()
                # Delete the leading and trailing blank characters in the line
                line = line.strip() 
                # Delete tab characetrs
                line = line.strip("\t")
                # Skip if the line contains nothing
                if (line == ""): continue
                # Divide the line into header, body
                header = ""
                body   = ""
                if (is_in_cmt == 1):
                    header, body = self.__get_header_cmt(line)
                    if (header):
                        is_in_cmt = 0
                else:
                    body = line
                # Check if there is an isolated, left open comment (*/)
                # [e.g.] aaaa /* bbbb */ cccc */ dddd
                # In this example, the last `*/` is syntax error.
                pos = self.__get_pos_left_open_cmt(body)
                if (pos >= 0):
                    msg = "An unclosed comment is detected!\n" \
                        + "Please check the following part of your code:" 
                    self.__print_error(msg)
                    print("{}".format(line))
                    sys.exit()
                # Replace a // comment by a /**/ comment if it exists
                body = self.__change_cmt_type(body) 
                # Delete /**/ comments that do not contains FDPS directives
                body = self.__delete_unnecessary_cmt(body)
                # Extract a footer comment if it exists
                body, footer = self.__get_footer_cmt(body)
                if (footer):
                    is_in_cmt = 1
                # Construct the reduced line
                line = header + body + footer
                if (line == ""):
                    continue
                # Update line_buf
                line_buf += line
                # Print for debugging
                #print("[dbg2] {0}".format(line))
                #print("header = {0}".format(header))
                #print("body   = {0}".format(body))
                #print("footer = {0}".format(footer))
                #print("line   = {0}".format(line))
                #print("line_buf = {0}".format(line_buf))
                #print("is_in_cmt = {0}".format(is_in_cmt))
                # Save the line
                if (is_in_cmt == 0):
                    line_buf = self.__delete_unnecessary_cmt(line_buf)
                    if (line_buf):
                        f.content.append(line)
                    line_buf = ""
            # Error handling
            if (is_in_cmt):
                # In this case, we rearch the end of the file although some /**/ comment
                # is not closed.
                self.__print_error("Some /**/ comment is not closed. Please check your code.")
                sys.exit()
            # Close the file
            fp.close()

    def analyze(self):
        # Output to stdout
        self.__print_checkpoint("analyzing C files ...")
        # Loop w.r.t. files
        num_left_brace = 0
        num_right_brace = 0
        for f in self.files:
            # Check the file name
            print("processing ... the file {0}".format(f.name))
            # Initialize flags and stacks
            struct_detect_mode = False
            num_left_curly_brace = 0
            num_right_curly_brace = 0
            #=============================================================================
            # Sweep the file content to detect structures
            print("identifying C structures ...")
            line_num = 0
            for line in f.content:
                #print("[dbg] {0}".format(line))
                len_line = len(line)
                start_srch = 0
                flag_re_srch = False # flag to control the while-loop below.
                while start_srch <= len_line - 1:
                    # [Note] This while-loop is required to deal with the case that
                    #        multiple structures are defined in a single line.
                    #print("struct_detect_mode = {0}".format(struct_detect_mode))
                    #print("start_srch = {0}".format(start_srch))
                    # Seek keyword `struct`
                    if (struct_detect_mode == False):
                        items = self.__find_struct_keyword(line, \
                                                           num_left_brace, \
                                                           num_right_brace, \
                                                           start_srch)
                        start_struct_key = items[0]
                        end_struct_key   = items[1]
                        num_left_brace   = items[2]
                        num_right_brace  = items[3] 
                        if (start_struct_key != -1):
                            # Create a struct object
                            f.structures.append(Structure_Data())
                            s = f.structures[-1]
                            # Setup struct
                            s.line_num_struct_key = line_num
                            s.start_struct_key = start_struct_key
                            s.end_struct_key = end_struct_key
                            # Update start_srch
                            start_srch = end_struct_key + 1
                            #print("start_srch       = {0}".format(start_srch))
                            #print("start_struct_key = {0}".format(start_struct_key))
                            #print("end_struct_key   = {0}".format(end_struct_key))
                            # Raise the flag
                            struct_detect_mode = True
                        else:
                            # In this case, we do not need to re-search the remaining text.
                            flag_re_srch = False

                    if (struct_detect_mode == True):
                        if (start_srch < len_line):
                            pos_non_cmt = self.__find_pos_non_cmt(line,start_srch)
                            if pos_non_cmt:
                                flag_break = False
                                while pos_non_cmt:
                                    start, end = pos_non_cmt.pop()
                                    for i in range(start,end+1):
                                        char = line[i]
                                        if (char == "{"):
                                            num_left_curly_brace += 1
                                            if (num_left_curly_brace == 1):
                                                s = f.structures[-1]
                                                s.line_num_left_curly_brace = line_num
                                                s.start_left_curly_brace = i
                                        elif (char == "}"):
                                            num_right_curly_brace += 1
                                            if (num_left_curly_brace > 0):
                                                if (num_left_curly_brace == num_right_curly_brace):
                                                    s = f.structures[-1]
                                                    s.line_num_right_curly_brace = line_num
                                                    s.start_right_curly_brace = i
                                                    # Reset 
                                                    num_left_curly_brace = 0
                                                    num_right_curly_brace = 0
                                                    # Update start_srch
                                                    start_srch = i + 1
                                                    #print("line     = {0}".format(line))
                                                    #print("len_line = {0}".format(len_line))
                                                    # Set the flags
                                                    flag_re_srch = True
                                                    flag_break   = True
                                                    struct_detect_mode = False
                                                    # Print for debugging
                                                    #print("struct is detected!")
                                                    break # escape from the for-loop
                                        # Check error
                                        if (num_right_curly_brace > num_left_curly_brace):
                                            msg = "Unparsable structure is detected!\n"
                                            self.__print_error(msg)
                                            sys.exit()
                                    if (flag_break):
                                        break # escape from the inner while-loop

                    # Check if we need to seek struct keyword again for this line
                    if (flag_re_srch == False):
                        break
                # Update line_num
                line_num += 1

            #=============================================================================
            # Check if the identified C structures are user-defined types or not
            print("checking if the identified C structure are eligible to be an user-defined type or not ...")
            line_num_lower_limit = 0
            for s in f.structures:
                # Print for debugging
                #print("--------")
                #for line_num in range(s.line_num_struct_key,s.line_num_right_curly_brace+1):
                #    print("[input] {0}".format(f.content[line_num]))
                #print("--------")
                # Reset flags
                is_eligible_to_be_fdps_UDT = False 
                # Check if this structure has tag name.
                # (note that tag name must be in between `struct` and '{'.)
                for line_num in range(s.line_num_struct_key,s.line_num_left_curly_brace+1):
                    line = f.content[line_num]
                    start_srch = 0
                    if (line_num == s.line_num_struct_key):
                        start_srch = s.end_struct_key + 1
                    end_srch = len(line) - 1
                    if (line_num == s.line_num_left_curly_brace):
                        end_srch = s.start_left_curly_brace - 1
                    start, end = self.__find_tag_name(line[:end_srch+1],start_srch) 
                    if (start != -1):
                        s.line_num_tag_name = line_num
                        s.start_tag_name = start
                        s.end_tag_name   = end
                        s.tag_name = line[start:end+1]
                        is_eligible_to_be_fdps_UDT = True
                        break
                # Check if there is a FDPS directive that specifies that
                # this C-structure corresponds to which user-defined type.
                if (is_eligible_to_be_fdps_UDT):
                    # First, check if there is a FDPS directive before `struct`.
                    line_num = s.line_num_struct_key
                    while line_num >= line_num_lower_limit:
                        line = f.content[line_num]
                        end_srch = len(line) - 1
                        if (line_num == s.line_num_struct_key):
                            end_srch = s.start_struct_key - 1
                        pos_cmt = self.__find_pos_cmt2(line[:end_srch+1])
                        if pos_cmt:
                            start, end = pos_cmt[-1] # Check the last comment only.
                            start += 2
                            end -= 2
                            s.fdps_str_dir = self.__get_fdps_str_dir(line[start:end+1])
                            break
                        line_num -= 1
                    # Next, check if there is a FDPS directive after tag name.
                    flag_break = False
                    for line_num in range(s.line_num_tag_name,s.line_num_right_curly_brace+1):
                        line = f.content[line_num]
                        start_srch = 0
                        if (line_num == s.line_num_tag_name):
                            start_srch = s.end_tag_name + 1
                        pos_cmt = self.__find_pos_cmt2(line)
                        if pos_cmt:
                            while pos_cmt:
                                start, end = pos_cmt.pop()
                                if (start >= start_srch):
                                    start += 2
                                    end -= 2
                                    if (s.fdps_str_dir):
                                        # In this case, this C-structure already has a FDPS directive.
                                        # Here, we check if there are multiple definitions or not.
                                        ret = self.__is_fdps_str_dir(line[start:end+1])
                                        if (ret):
                                            msg = "It is not allowed to describe multiple FDPS directives that\n" \
                                                + "specifies that this C-structre corresponds to which user-defined type!\n" 
                                            self.__print_error(msg)
                                            msg = "The script were processing the structure of the name of {0}".format(s.tag_name)
                                            print(msg)
                                            sys.exit()
                                    else:
                                        s.fdps_str_dir = self.__get_fdps_str_dir(line[start:end+1])
                                    flag_break = True
                                    break # break the inner while loop.
                        if (flag_break):
                            break # break the outer for loop.

                # Check if the current C-structure is eligible to be user-defined type or not.
                # If not, skip the remaining part of the for-loop.
                line_num_lower_limit = s.line_num_right_curly_brace
                if (not s.fdps_str_dir):
                    continue

            #=============================================================================
            print("deleting the information about C structures which are not user-defined types ...")
            structures_tmp = f.structures
            f.structures = []
            for s in structures_tmp:
                if (s.fdps_str_dir):
                    f.structures.append(s)
                else:
                    if (s.tag_name):
                        msg = "struct {0} is not eligible to be an user-defined type. ".format(s.tag_name) \
                            + "Hence, deleted!"
                    else:
                        msg = "strcutre without having a tag name cannot be an user-defined type. " \
                            + "Hence, deleted!"
                    self.__print_warning(msg)
            structures_tmp = [] 

            #=============================================================================
            # Check if the identified C structures are user-defined types or not
            print("identify other FDPS directives and member variables from user-defined types ...")
            for s in f.structures:
                # Identify candidates of FDPS directive for specifying the ways of data processing
                # in FDPS. (note that full check of the format and syntax is performed later)
                for line_num in range(s.line_num_left_curly_brace,s.line_num_right_curly_brace+1):
                    line = f.content[line_num]
                    start_srch = 0
                    if (line_num == s.line_num_left_curly_brace):
                        start_srch = s.start_left_curly_brace + 1
                    pos_cmt = self.__find_pos_cmt2(line)
                    if pos_cmt:
                        while pos_cmt:
                            start, end = pos_cmt.pop()
                            if (start >= start_srch):
                                start += 2
                                end -= 2
                                ret = self.__get_fdps_meth_dir_cand(line[start:end+1])
                                if (ret):
                                    s.fdps_meth_dirs.append(ret)

                # Print for debugging 
                #for line_num in range(s.line_num_struct_key,s.line_num_right_curly_brace+1):
                #    line = f.content[line_num]
                #    print("[identified.] {0}".format(line))
                                        
                # Identify and extract the variable declaration parts
                num_left_curly_brace     = 0
                num_right_curly_brace    = 0
                is_1st_normal_char_found = False
                is_end_of_var_decl_part  = False
                line_num_1st_normal_char = -1
                pos_1st_normal_char      = -1
                line_num_semicolon       = -1
                pos_semicolon            = -1
                text                     = ""
                for line_num in range(s.line_num_left_curly_brace,s.line_num_right_curly_brace+1):
                    line = f.content[line_num]
                    start_srch = 0
                    if (line_num == s.line_num_left_curly_brace):
                        start_srch = s.start_left_curly_brace + 1
                        if (start_srch > len(line) - 1):
                            continue
                    end_srch = len(line) - 1
                    if (line_num == s.line_num_right_curly_brace):
                        end_srch = s.start_right_curly_brace - 1
                        if (end_srch < 0):
                            continue
                    # Read the non-comments parts of the text.
                    pos_non_cmt = self.__find_pos_non_cmt(line,start_srch,end_srch)
                    if pos_non_cmt:
                        while pos_non_cmt:
                            start, end = pos_non_cmt.pop()
                            for i in range(start,end+1):
                                char = line[i]
                                if (char == "{"):
                                    num_left_curly_brace += 1
                                elif (char == "}"):
                                    num_right_curly_brace += 1
                                elif (char == ";"):
                                    if (num_left_curly_brace == num_right_curly_brace):
                                        line_num_semicolon      = line_num
                                        pos_semicolon           = i
                                        is_end_of_var_decl_part = True
                                if (is_1st_normal_char_found == False):
                                    m = re.match(r"[a-zA-Z_]",char)
                                    if (m):
                                        line_num_1st_normal_char = line_num
                                        pos_1st_normal_char      = i
                                        is_1st_normal_char_found = True
                                text += char
                                if (is_end_of_var_decl_part):
                                    # Create a var. decl. object
                                    s.var_decl_parts.append(Variable_Declaration_Data())
                                    vdp = s.var_decl_parts[-1]
                                    # Set var. decl. object
                                    vdp.line_num_start = line_num_1st_normal_char 
                                    vdp.pos_start      = pos_1st_normal_char
                                    vdp.line_num_end   = line_num_semicolon
                                    vdp.pos_end        = pos_semicolon
                                    vdp.text           = text
                                    # Reset auxiliary variables
                                    num_left_curly_brace     = 0
                                    num_right_curly_brace    = 0
                                    is_1st_normal_char_found = False
                                    is_end_of_var_decl_part  = False
                                    line_num_1st_normal_char = -1
                                    pos_1st_normal_char      = -1
                                    line_num_semicolon       = -1
                                    pos_semicolon            = -1
                                    text                     = ""
                # Print for debugging
                #for k in range(0,len(s.var_decl_parts)):
                #    vdp = s.var_decl_parts[k]
                #    char_1st_normal = f.content[vdp.line_num_start][vdp.pos_start]
                #    char_semicolon  = f.content[vdp.line_num_end][vdp.pos_end]
                #    print("[chk] {0} {1} {2}".format(char_1st_normal,char_semicolon,vdp.text))

                # Identify member variables
                for vdp in s.var_decl_parts:
                    num_left_curly_brace  = 0
                    num_right_curly_brace = 0
                    is_able_to_split      = False
                    is_1st_split          = True
                    data_type             = ""
                    var_names             = []
                    for char in vdp.text:
                        if (char == "{"):
                            num_left_curly_brace += 1
                        elif (char == "}"):
                            num_right_curly_brace += 1
                        elif ((char == ",") or (char == ";")):
                            if (num_left_curly_brace == num_right_curly_brace):
                                is_able_to_split = True
                        if (is_able_to_split):
                            if (is_1st_split):
                                # In this case, we extract data type.
                                items = text.split() # split by blank characters
                                if (len(items) < 2):
                                    msg = "The name of either data type or variable name is missing!"
                                    self.__print_error(msg)
                                    sys.exit()
                                # Extract the name of data type.
                                # (here we assume that the last item represents the name of variable.)
                                data_type = ""
                                for i in range(0,len(items)-1):
                                    data_type += items[i].strip() + " "
                                data_type = data_type.strip()
                                # Check if the name of data type is valid or not
                                is_valid_data_type = True
                                for key in data_type.split():
                                    if key not in self.__usable_data_types:
                                        is_valid_data_type = False
                                if (is_valid_data_type == False):
                                    msg = "Data type `{0}` cannot be used for a FDPS user-defined type.".format(data_type) 
                                    self.__print_error(msg)
                                    sys.exit()
                                # Get the variable name
                                var_name = items[-1].strip()
                                if var_name in self.__reserved_words_in_c:
                                    msg = "Cannot use a reserved word in C for the name of a variable!"
                                    self.__print_error(msg)
                                    sys.exit()
                                var_names.append(var_name)
                                # Update is_1st_split 
                                is_1st_split = False
                            else:
                                # Get the variable name
                                var_names.append(text.strip())
                            # Update auxiliary variables
                            text             = ""
                            is_able_to_split = False
                        else:
                            text += char
                    # Print for debugging
                    #print("[chk] ",data_type,var_names)

                    # Check if there are FDPS directives for this (these) member variables
                    # (note that only one FDPS directive is allowed.)
                    attrib = ""
                    # (1) Check if there is a FDPS directive in the previous line
                    line = f.content[vdp.line_num_start-1]
                    pos_cmt = self.__find_pos_cmt2(line)
                    if pos_cmt:
                        if len(pos_cmt) > 1:
                            msg = "Only one FDPS directive is allowed in a line!"
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            start, end = pos_cmt.pop()
                            if (start != 0):
                                # In this case, there is no isolated FDPS directive
                                # in the previous line. There is also a variable
                                # declaration part in the previous line.
                                pass
                            else:
                                start += 2
                                end -= 2
                                attrib = self.__get_fdps_mbr_dir(line[start:end+1])
                    else:
                        # In this case, there are no FDPS directive in the previous line
                        pass
                    # (2) Check if there is a FDPS directive after the character `;`.
                    line = f.content[vdp.line_num_end]
                    if (vdp.pos_end + 1 < len(line) - 1):
                        pos_cmt = self.__find_pos_cmt2(line[vdp.pos_end+1:])
                        if pos_cmt:
                            if len(pos_cmt) > 1:
                                msg = "Only one FDPS directive is allowed in a line!"
                                self.__print_error(msg)
                                sys.exit()
                            else:
                                start, end = pos_cmt.pop()
                                start += (vdp.pos_end + 1) + 2
                                end   += (vdp.pos_end + 1) - 2
                                attrib = self.__get_fdps_mbr_dir(line[start:end+1])
                    # Create Member_Data objects 
                    for var_name in var_names:
                        # Check if var_name is array or not
                        is_array  = False
                        array_dim = -1
                        m = re.match(r"^[a-zA-Z_].*",var_name)
                        if (m):
                            # In this case, the name of the variable starts with
                            # a valid character. Next, we check this variable is
                            # a scalar or an array.
                            m = re.match(r"^[a-zA-Z_][a-zA-Z0-9_]+$",var_name)
                            if (m):
                                # In this case, this variable must be a scalar.
                                pass
                            else:
                                # In this case, this variable may be an array.
                                # In order to be an array, the following conditions 
                                # must be fulfilled:
                                # (1) the name of variable must consist of the set
                                #     of characters [a-zA-Z0-9_[\]].
                                # (2) `[` and `]` must be grouped into pairs.
                                # (3) strings of characters in the brackets must be 
                                #     integral numbers.
                                # In addition, by the restriction from Fortran i/f,
                                # (4) it must be an 1D array.
                                m = re.match(r"^[a-zA-Z0-9_[\]]*$",var_name)
                                if (m):
                                    # Condition (1) is fulfilled in this case.
                                    # Next, we check if Condition (2) is fulfilled.
                                    num_left_bracket = 0
                                    num_right_bracket = 0
                                    for char in var_name:
                                        if (char == "["):
                                            num_left_bracket += 1
                                        elif (char == "]"):
                                            num_right_bracket += 1
                                    if (num_left_bracket != num_right_bracket):
                                        msg = "There are unpaired brackets!" 
                                        self.__print_error(msg)
                                        sys.exit()
                                    # Then, we check if Condition (3) is fulfilled.
                                    buf = var_name
                                    while buf:
                                        m = re.search(r"[[](?P<array_dim>.*?)[\]]",buf)
                                        if (m):
                                            array_dim = m.group('array_dim')
                                            buf = buf[m.end():]
                                            m = re.match(r"[0-9]+",array_dim)
                                            if (m):
                                                # In this case, the array size are described by
                                                # integral numbers.
                                                pass
                                            else:
                                                msg = "Non-integral numbers are used to specify the size of array!"
                                                self.__print_error(msg)
                                                sys.exit()
                                        else:
                                            break
                                    is_array = True
                                    # Finally, we check if Condition (4) is fulfilled.
                                    if (num_left_bracket != 1):
                                        msg = "Currently, only 1D array is allowed in a FDPS user-defined type!"
                                        self.__print_error(msg)
                                        sys.exit()
                                    # Extract true variable name
                                    m = re.match(r"[a-zA-Z_][a-zA-Z0-9_]+",var_name)
                                    var_name = m.group()
                                else:
                                    msg = "Invalid characters are used for variable name!"
                                    self.__print_error(msg)
                                    sys.exit()

                        else:
                            msg = "The name of member variable starts with an invalid character!"
                            self.__print_error(msg)
                            sys.exit()
                        # Registier data
                        s.members.append(Member_Data())
                        mbr = s.members[-1]
                        mbr.name      = var_name
                        mbr.attrib    = attrib
                        mbr.data_type = data_type
                        mbr.is_array  = is_array
                        mbr.array_dim = array_dim
                # Print for debugging
                #for mbr in s.members:
                #    print("{0}, {1}, {2}, {3}, {4}".format(mbr.data_type,mbr.name,mbr.is_array,mbr.array_dim,mbr.attrib))

    def check(self):
        # Output to stdout
        self.__print_checkpoint("checking if each structure fulfills our requirements...")
        # Check if each structure fulfills requirements of FDPS:
        # (1) Check if all characters in tag name are lower case.
        # (2) Check if name collision occurs if all characters in the names of member
        #     variables are decapitalized..
        # (3) Check if member variables more than two have same FDPS directive 
        #     like `$fdps position` ?
        # (4) Check if data types of member variables having a FDPS directive are
        #     consistent with those directives.
        for f in self.files:
            for s in f.structures:
                # Perform checking (1)
                rslt = []
                for char in s.tag_name:
                    if char != '_':
                        rslt.append(char.islower())
                if not all(rslt):
                    msg = """
                    tag name of strcture {0} contains capital letters!

                    Information
                    ===========
                       In the current version of the C i/f to FDPS, all letters in tag
                       name must be lower-case in order for this structure to be FDPS's
                       user-defined type.
                    """.format(s.tag_name)
                    msg = msg[1:].rstrip() + "\n"
                    msg = textwrap.dedent(msg)
                    self.__print_error(msg)
                    sys.exit()
                # Perfrom checking (2)
                mbr_set = set()
                for mbr in s.members:
                    mbr.name = mbr.name.lower() 
                    mbr_set.add(mbr.name)
                if (len(s.members) != len(mbr_set)):
                    msg = """
                    Names of member variables are duplicated in structure {0}!

                    Information
                    ===========
                       In the current version of the C i/f to FDPS, all of member variable
                       names are automatically decapitalized. After this decapitalization,
                       structure must fulfill the grammer of the C language.
                    """.format(s.tag_name)
                    msg = msg[1:].rstrip() + "\n"
                    msg = textwrap.dedent(msg)
                    self.__print_error(msg)
                    sys.exit()
                # Perform checking (3)
                attrib_list = []
                attrib_set  = set()
                for mbr in s.members:
                    if (mbr.attrib != ""):
                        attrib_list.append(mbr.attrib)
                        attrib_set.add(mbr.attrib)
                if (len(attrib_list) != len(attrib_set)):
                    msg = """
                    There are duplicative FDPS directives in structure {0}!
                    """.format(s.tag_name)
                    msg = msg[1:].rstrip() + "\n"
                    msg = textwrap.dedent(msg)
                    self.__print_error(msg)
                    sys.exit()
                # Perform checking (4)
                for mbr in s.members:
                    if (mbr.attrib == "id"):
                        if not self.__is_s64(mbr):
                            msg = """
                            The data type of member variable {0} of strcture {1} must be type fdps_s64!
                            """.format(mbr.name, s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (mbr.attrib == "charge" or mbr.attrib == "rsearch"):
                        if not (self.__is_f32(mbr) or self.__is_f64(mbr)):
                            msg = """
                            The data type of member variable {0} of structure {1} must be floating
                            point number types!
                            """.format(mbr.name, s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (mbr.attrib == "position" or mbr.attrib == "velocity"):
                        if not (self.__is_f32vec(mbr) or self.__is_f64vec(mbr)):
                            msg = """
                            The data type of member variable {0} of structure {1} must be vector type!
                            """.format(mbr.name, s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                # Register
                m = re.match(r"\$fdps\s+(?P<keywords>\S+)",s.fdps_str_dir)
                keys = m.group('keywords').split(',')
                for key in keys:
                    key = key.strip()
                    if (key == "FP"):
                        if not s.tag_name in self.__FPs:
                            self.__FPs.append(s.tag_name)
                        else:
                            msg = """
                            Structure {0} cannot be registered as FP because a structure 
                            having the same tag name is already defined.
                            """.format(s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (key == "EPI"):
                        if not s.tag_name in self.__EPIs:
                            self.__EPIs.append(s.tag_name)
                        else:
                            msg = """
                            Structure {0} cannot be registered as EPI because a structure 
                            having the same tag name is already defined.
                            """.format(s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (key == "EPJ"):
                        if not s.tag_name in self.__EPJs:
                            self.__EPJs.append(s.tag_name)
                        else:
                            msg = """
                            Structure {0} cannot be registered as EPJ because a structure 
                            having the same tag name is already defined.
                            """.format(s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (key == "Force"):
                        if not s.tag_name in self.__Forces:
                            self.__Forces.append(s.tag_name)
                        else:
                            msg = """
                            Structure {0} cannot be registered as EPJ because a structure 
                            having the same tag name is already defined.
                            """.format(s.tag_name)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()

        self.__print_checkpoint("performing syntax checking of FDPS directives ...")
        # Check the format and syntax of FDPS directives for specifying the ways of data
        # processing in FDPS.
        # For copyFromFP, copyFromForce, copyFromForcePM:
        # (1) Check if structure described in directive does exist?
        # (2) Check if names of member variables described in directive do exist?
        # (3, TODO) Check if the data types of the members in each pair are consistent
        #           with each other? (only for copyFromFP & copyFromForce)
        # For clear
        # (1) Check if directive described by a user fulfills the syntax?
        # (2) Convert initial values described in C manner into in Fortran manner.
        for f in self.files:
            for s in f.structures:
                for meth_dir_num in range(0,len(s.fdps_meth_dirs)):
                    meth_dir = s.fdps_meth_dirs[meth_dir_num]
                    items = meth_dir.split()
                    meth_type = items[1]
                    if (meth_type == "copyFromFP"):
                        tag_name_src = items[2]
                        if not self.__struct_existence_check(tag_name_src):
                            msg = """
                            Structure of the name of `{0}` that is specified as the copy source
                            in the following FDPS directive does not exist:
                            {1}
                            """.format(tag_name_src, meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                        if not tag_name_src in self.__FPs:
                            msg = """
                            Structure of the name of `{0}` is not FullParticle type!
                            Please check the description of the directive:
                            {1}
                            """.format(tag_name_src, meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                        text = ""
                        for i in range(3,len(items)):
                            text += items[i]
                        text = text.strip()
                        start_srch = 0
                        while start_srch <= len(text) - 1:
                            m = re.match(r"(?P<pair_info>[(].*?[)])",text[start_srch:])
                            if (m):
                                pair_info = m.group('pair_info')
                                start_srch += m.end('pair_info')
                                m = re.match(r"[(](?P<src>[a-zA-Z_][a-zA-Z0-9_]+),(?P<dest>[a-zA-Z_][a-zA-Z0-9_]+)[)]",pair_info)
                                if (m):
                                    mbr_name_src = m.group('src')
                                    mbr_name_dest = m.group('dest')
                                    if not self.__member_existence_check(tag_name_src,mbr_name_src):
                                        msg = """
                                        Structure `{0}` does not have a member variable of the name of `{1}`!
                                        Please check the description of the directive:
                                        {2}
                                        """.format(tag_name_src, mbr_name_src, meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                    if not self.__member_existence_check(s.tag_name,mbr_name_dest):
                                        msg = """
                                        Structure `{0}` does not have a member variable of the name of `{1}`!
                                        Please check the description of the directive:
                                        {2}
                                        """.format(s.tag_name, mbr_name_dest, meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                else:
                                    msg = """
                                    Syntax error is detected in the following FDPS directive:
                                    {0}
                                    """.format(meth_dir)
                                    msg = msg[1:].rstrip() + "\n"
                                    msg = textwrap.dedent(msg)
                                    self.__print_error(msg)
                                    sys.exit()
                    elif (meth_type == "copyFromForce"):
                        tag_name_src = items[2]
                        if not self.__struct_existence_check(tag_name_src):
                            msg = """
                            Structure of the name of `{0}` that is specified as the copy source
                            in the following FDPS directive does not exist:
                            {1}
                            """.format(tag_name_src, meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                        if not tag_name_src in self.__Forces:
                            msg = """
                            Structure of the name of `{0}` is not FORCE type!
                            Please check the description of the directive:
                            {1}
                            """.format(tag_name_src, meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                        text = ""
                        for i in range(3,len(items)):
                            text += items[i]
                        text = text.strip()
                        start_srch = 0
                        while start_srch <= len(text) - 1:
                            m = re.match(r"(?P<pair_info>[(].*?[)])",text[start_srch:])
                            if (m):
                                pair_info = m.group('pair_info')
                                start_srch += m.end('pair_info')
                                m = re.match(r"[(](?P<src>[a-zA-Z_][a-zA-Z0-9_]+),(?P<dest>[a-zA-Z_][a-zA-Z0-9_]+)[)]",pair_info)
                                if (m):
                                    mbr_name_src = m.group('src')
                                    mbr_name_dest = m.group('dest')
                                    if not self.__member_existence_check(tag_name_src,mbr_name_src):
                                        msg = """
                                        Structure `{0}` does not have a member variable of the name of `{1}`!
                                        Please check the description of the directive:
                                        {2}
                                        """.format(tag_name_src, mbr_name_src, meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                    if not self.__member_existence_check(s.tag_name,mbr_name_dest):
                                        msg = """
                                        Structure `{0}` does not have a member variable of the name of `{1}`!
                                        Please check the description of the directive:
                                        {2}
                                        """.format(s.tag_name, mbr_name_dest, meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                else:
                                    msg = """
                                    Syntax error is detected in the following FDPS directive:
                                    {0}
                                    """.format(meth_dir)
                                    msg = msg[1:].rstrip() + "\n"
                                    msg = textwrap.dedent(msg)
                                    self.__print_error(msg)
                                    sys.exit()
                    elif (meth_type == "copyFromForcePM"):
                        mbr_name = items[2]
                        if not self.__member_existence_check(s.tag_name, mbr_name):
                            msg = """
                            Structure `{0}` does not have a member variable of the name of `{1}`!
                            Please check the description of the directive:
                            {2}
                            """.format(s.tag_name, mbr_name, meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                        if (len(items) > 3):
                            msg = """
                            Variable names more than two are specified in copyFromForcePM directive!
                            Please check the description of the directive:
                            {0}
                            """.format(meth_dir)
                            msg = msg[1:].rstrip() + "\n"
                            msg = textwrap.dedent(msg)
                            self.__print_error(msg)
                            sys.exit()
                    elif (meth_type == "clear"):
                        if (len(items) > 2):
                            # In this case, the following two cases hold.
                            # $fdps clear mbr=val, ...
                            # $fdps clear function function_name
                            if (items[2] == "function"):
                                # In this case, an external function is used to initialize FORCE type.
                                # Here, we check if its name fulfills the C grammer.
                                meth_dir_new = "$fdps clear subroutine "
                                if len(items) != 4:
                                    msg = """
                                    Syntax error is detected in clear directive:
                                    {0}
                                    The directive must be written as `$fdps clear function function_name`.
                                    """.format(meth_dir)
                                    msg = msg[1:].rstrip() + "\n"
                                    msg = textwrap.dedent(msg)
                                    self.__print_error(msg)
                                    sys.exit()
                                func_name = items[3]
                                m = re.match(r"^[a-zA-Z_][a-zA-Z0-9_]+$",func_name)
                                if (m):
                                    # Because of the specification of Fortran i/f, all the letters
                                    # in the function name must be lower case. So, we check this here.
                                    m = re.match(r"^[a-z_][a-z0-9_]+$",func_name)
                                    if (m):
                                        # In this case, all the requirements are fulfilled.
                                        # We can make the Fortran version of the directive.
                                        meth_dir_new += func_name
                                        s.fdps_meth_dirs[meth_dir_num] = meth_dir_new
                                    else:
                                        msg = """
                                        A capital letter is found in the function name in the following
                                        clear directive:
                                        {0}

                                        Information
                                        ===========
                                           In the current version of the C i/f to FDPS, all of letters
                                           in function name of clear directive must be lower case.
                                        """.format(meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                else:
                                    msg = """
                                    A function name that is specified in the following clear directive
                                    is in violation with the C grammer:
                                    {0}
                                    """.format(meth_dir)
                                    msg = msg[1:].rstrip() + "\n"
                                    msg = textwrap.dedent(msg)
                                    self.__print_error(msg)
                                    sys.exit()
                            else:
                                # In this case, we replace initial values described in C manner
                                # to those described in Fortran manner.
                                meth_dir_new = "$fdps clear "
                                # First, get a part of the directive after "$fdps clear".
                                m = re.search(r"\$fdps\s*clear",meth_dir)
                                text = meth_dir[m.end():].strip()
                                # Then, split the obtained text into each initialization instruction.
                                # When spliting the text, we have to treat the problems.
                                # (i) A scalar variable or array of char-type may be initialized by 
                                #     a string constant that contains a commna (string constant is defined
                                #     as characters surrounded by double quotes). Therefore, we cannot
                                #     use Python's split method.
                                # (ii) Users may try to initialize char-type variables by a CHARACTER
                                #      CONSTANT (defined as a single character surrounded by single
                                #      quotes). In the C++ side, we initialize char-type variable by
                                #      strcpy function, which does not accept a character constant
                                #      (Segmentation fault occurs in run-time).
                                # (iii) Also, uses may try to initialize char-type variables by
                                #       characters surrounded by single quotes, which violates the 
                                #       C grammer.
                                # This script deals with the problems by the following way.
                                # Read the text from the head character by character with counting
                                # the numbers of single and double quotes. If we find a commna, we check
                                # if the numbers are even numbers or not. If those are even, we can
                                # adopt the comann as delimiter.
                                # If we detect cases (ii) and (iii), we automatically convert them into
                                # string constats.
                                num_sng_quote = 0
                                num_dbl_quote = 0
                                init_insts = []
                                i_prev = -1
                                for i in range(0,len(text)):
                                    char = text[i]
                                    if (char == "\""):
                                        if (num_sng_quote % 2 == 0):
                                            num_dbl_quote += 1
                                    elif (char == "'"):
                                        if (num_dbl_quote % 2 == 0):
                                            num_sng_quote += 1
                                    elif (char == ","):
                                        if ((num_sng_quote % 2 == 0) and \
                                            (num_dbl_quote % 2 == 0)):
                                            init_insts.append(text[i_prev+1:i].strip())
                                            i_prev = i
                                if ((num_sng_quote % 2 == 0) and \
                                    (num_dbl_quote % 2 == 0)):
                                    init_insts.append(text[i_prev+1:].strip())
                                else:
                                    msg = """
                                    the following clear directive contains unclosed pairs of single
                                    or double quotes. Please check the directive:
                                    {0}
                                    """.format(meth_dir)
                                    msg = msg[1:].rstrip() + "\n"
                                    msg = textwrap.dedent(msg)
                                    self.__print_error(msg)
                                    sys.exit()
                                init_insts.reverse() 
                                # This reverse is to keep the original order after taking the items by pop.
                                # Print for debugging --- [start] ---
                                #print("num_sng_quote = {0}".format(num_sng_quote))
                                #print("num_dbl_quote = {0}".format(num_dbl_quote))
                                #print(init_insts)
                                #sys.exit()
                                # --- [end] ----
                                is_1st=True
                                while init_insts:
                                    init_inst = init_insts.pop().strip()
                                    m = re.match(r"^(?P<mbr_name>[a-zA-Z_][a-zA-Z0-0_]+)\s*=(?P<val>.+)$",init_inst)
                                    if (m):
                                        mbr_name = m.group('mbr_name').strip()
                                        val = m.group('val').strip()
                                        if not self.__member_existence_check(s.tag_name,mbr_name):
                                            msg = """
                                            Structure {0} does not have a member variable of the name of `{0}`!
                                            Please check the directive:
                                            {2}
                                            """.format(s.tag_name, mbr_name, meth_dir)
                                            msg = msg[1:].rstrip() + "\n"
                                            msg = textwrap.dedent(msg)
                                            self.__print_error(msg)
                                            sys.exit()
                                        if val == "keep":
                                            if not is_1st:
                                                meth_dir_new += ", "
                                            meth_dir_new += "{0}=keep".format(mbr_name)
                                            is_1st = False
                                            continue
                                        data_type = self.__get_data_type(s.tag_name, mbr_name)
                                        if self.__is_integer(data_type):
                                            m = re.match(r"^(?P<num>(\+\s*|-\s*)?[0-9]+)([uUlL]|[uU][lL]|[lL][uU])?$",val)
                                            if (m):
                                                val = m.group('num')
                                                if not is_1st:
                                                    meth_dir_new += ", "
                                                meth_dir_new += "{0}={1}".format(mbr_name, val)
                                                is_1st = False
                                            else:
                                                msg = """
                                                Non-integer value `{0}` cannot be substituted into member variable `{1}`!
                                                Please check the directive:
                                                {2}

                                                Information
                                                ===========
                                                   In the current version of the C i/f to FDPS, a numerical expression such as
                                                   1+2, 1.0f-1 (this means 1.0f - 1.0 because f is the qualifier to represent
                                                   the accuracy of data type), 1.0/2.0, etc. cannot be used for an initial value.
                                                """.format(val, mbr_name, meth_dir)
                                                msg = msg[1:].rstrip() + "\n"
                                                msg = textwrap.dedent(msg)
                                                self.__print_error(msg)
                                                sys.exit()
                                        elif self.__is_float(data_type):
                                            m_float = re.match(r"^(\+\s*|-\s*)?[0-9]*\.[0-9]*(e-|e|e\+|E-|E|E\+)?[0-9]*[fFlL]?$",val)
                                            m_int   = re.match(r"^(\+\s*|-\s*)?[0-9]*$",val)
                                            if (m_float or m_int):
                                                if (m_float):
                                                    # In this case, we replace f or F by e, and e or E by d.
                                                    # In addition, we delete the suffix l or L if it exists.
                                                    val = re.sub(r"[eE]",r"d",val)
                                                    val = re.sub(r"[fF]",r"e",val)
                                                    val = re.sub(r"[lL]",r"",val)
                                                if not is_1st:
                                                    meth_dir_new += ", "
                                                meth_dir_new += "{0}={1}".format(mbr_name, val)
                                                is_1st = False
                                            else:
                                                msg = """
                                                Non-floating point value `{0}` cannot be substituted into member variable `{1}`!
                                                Please check the directive:
                                                {2}

                                                Information
                                                ===========
                                                   In the current version of the C i/f to FDPS, a numerical expression such as
                                                   1+2, 1.0f-1 (this means 1.0f - 1.0 because f is the qualifier to represent
                                                   the accuracy of data type), 1.0/2.0, etc. cannot be used for an initial value.
                                                """.format(val, mbr_name, meth_dir)
                                                msg = msg[1:].rstrip() + "\n"
                                                msg = textwrap.dedent(msg)
                                                self.__print_error(msg)
                                                sys.exit()
                                        elif self.__is_boolean(data_type):
                                            m = re.match(r"^true$|^false$",val)
                                            if (m):
                                                val = re.sub(r"true",r".true.",val)
                                                val = re.sub(r"false",r".false.",val)
                                                if not is_1st:
                                                    meth_dir_new += ", "
                                                meth_dir_new += "{0}={1}".format(mbr_name, val)
                                                is_1st = False
                                            else:
                                                msg = """
                                                Non-boolean value `{0}` cannot be substituted into member variable `{1}`!
                                                Please check the directive:
                                                {2}
                                                """.format(val, mbr_name, meth_dir)
                                                msg = msg[1:].rstrip() + "\n"
                                                msg = textwrap.dedent(msg)
                                                self.__print_error(msg)
                                                sys.exit()
                                        elif self.__is_character(data_type):
                                            m = re.match(r"^\".*\"$|^'.*'$",val)
                                            if (m):
                                                m = re.match(r"^'.*'$",val)
                                                if (m):
                                                    val = re.sub(r"^'","\"",val)
                                                    val = re.sub(r"'$","\"",val)
                                                if not is_1st:
                                                    meth_dir_new += ", "
                                                meth_dir_new += "{0}={1}".format(mbr_name, val)
                                                is_1st = False
                                            else:
                                                msg = """
                                                Non string of characters `{0}` cannot be substituted into member variable `{1}`!
                                                Please check the directive:
                                                {2}
                                                """.format(val, mbr_name, meth_dir)
                                                msg = msg[1:].rstrip() + "\n"
                                                msg = textwrap.dedent(msg)
                                                self.__print_error(msg)
                                                sys.exit()
                                    else:
                                        msg = """
                                        Syntax error is detected in the clear directive!
                                        clear directive must be written like this:
                                        $fdps clear mbr0=val0, mbr1=val1, ...
                                        Please check the directive:
                                        {0}
                                        """.format(meth_dir)
                                        msg = msg[1:].rstrip() + "\n"
                                        msg = textwrap.dedent(msg)
                                        self.__print_error(msg)
                                        sys.exit()
                                # Substitute a replaced clear directive
                                s.fdps_meth_dirs[meth_dir_num] = meth_dir_new
                                # Print for debugging --- [start] ---
                                #print("meth_dir (repl) = {0}".format(meth_dir_new))
                                #sys.exit()
                                # --- [end] ---


    def generate(self,output_dir):
        # Output to stdout
        self.__print_checkpoint("generating FDPS C interface programs ....")
        # Make the directory
        output_dir = output_dir.strip()
        if (os.path.exists(output_dir) is False):
            os.makedirs(output_dir)
        # Make the working directory
        work_dir = ""
        while 1:
            work_dir = "./" + self.__get_random_string(16)
            if not os.path.isdir(work_dir):
                os.makedirs(work_dir)
                print("A work directory {0} is created!".format(work_dir))
                break
        # Compute blueprint_dir 
        blueprint_dir = os.path.abspath(os.path.dirname(__file__)) \
                      + "/../src/c_interface/blueprints"
        # Generate interface programs
        #---------------------------------------------------------------------------------
        # (1) main.cpp
        input_file  = blueprint_dir + "/main_blueprint.cpp"
        output_file = output_dir    + "/main.cpp"
        ifh = open(input_file,'r')
        ofh = open(output_file,'w')
        for line in ifh:
            ofh.write(line) # simply copy
        ifh.close()
        ofh.close()
        #---------------------------------------------------------------------------------
        # (2) user_defined.F90
        file_name = work_dir + "/user_defined.F90"
        fh = open(file_name,"w")
        # Write the header of the files
        text = """
        !===============================
        !   MODULE: User defined types
        !===============================
        module user_defined_types
           use, intrinsic :: iso_c_binding
           use fdps_vector
           use fdps_super_particle
           implicit none
        """
        text = text[1:].rstrip() + "\n\n"
        text = textwrap.dedent(text)
        fh.write(text)
        # Write the definitions of Fortran derived data types
        for f in self.files:
            for s in f.structures:
                # Skip non-FDPS user defined types
                if (s.fdps_str_dir == ""):
                    continue
                # Write the header of the class
                fdps_dir = "!" + s.fdps_str_dir.strip()
                fh.write("   type, public, bind(c) :: {0} {1}\n".format(s.tag_name, \
                                                                        fdps_dir))
                # Write the FDPS directive for methods
                for c in s.fdps_meth_dirs:
                    fdps_dir = "      !" + c.strip() + "\n"
                    fh.write(fdps_dir)
                # Write the definitions of the member variables
                for mbr in s.members:
                    # Get data type in FDPS-style
                    data_type, array_dim = self.__get_fdpsDT(mbr)
                    # FDPS directive for member variables
                    fdps_dir = ""
                    if (mbr.attrib):
                        fdps_dir = "!$fdps " + mbr.attrib
                    if (array_dim != 0):
                        fh.write("      {0} :: {1}({2}) {3}\n".format(data_type, \
                                                                      mbr.name,  \
                                                                      array_dim, \
                                                                      fdps_dir))
                    else:
                        fh.write("      {0} :: {1} {2}\n".format(data_type, \
                                                                 mbr.name,  \
                                                                 fdps_dir))
                # Write the footer of the class
                fh.write("   end type {0}\n".format(s.tag_name))
                # Write black lines
                fh.write("\n\n")
        # Write the footer of the file
        fh.write("end module user_defined_types")
        # Close the file
        fh.close()
        #---------------------------------------------------------------------------------
        # Perform gen_ftn_if.py
        script_name    = "gen_ftn_if.py"
        script_options = "-fmode-to-gen-c-if"
        script_dir     = os.path.abspath(os.path.dirname(__file__))
        cmd = script_dir + "/" + script_name + " " \
            + work_dir + "/user_defined.F90" + " " \
            + script_options
        subprocess.call(cmd,shell=True)
        # Delete the working directory
        shutil.rmtree(work_dir)
        print("Working directory {0} is deleted!".format(work_dir))


#================================
#   Function definition
#================================
def analyze_CL_args():
    #-------------------------------------------------------------------------------------
    # More information about the module `argparser` can be found at the following sites:
    #  Python 3.5) http://docs.python.jp/3.5/library/argparse.html
    #  Python 2.x) http://docs.python.jp/2/library/argparse.html
    #-------------------------------------------------------------------------------------
    # Create an ArgumentParser object
    description = "Analyze user's C codes and generate C++/C source files " \
                + "required to use FDPS from the user's C code."
    parser = argparse.ArgumentParser(description=description)
    # Input option
    help_msg = "The PATHs of input C files"
    parser.add_argument('input_files', nargs="+", \
                        type=str, \
                        help=help_msg, metavar='FILE')
    # Output option
    help_msg = "The PATH of output directory"
    parser.add_argument('-o','--output','--output_dir', \
                        default="./",type=str,required=False, \
                        help=help_msg, metavar='DIRECTORY', \
                        dest='output_dir')
    # Space dimension option 
    help_msg = "Indicate that simulation is performed in the 2-dimensional space " \
             + "(equivalent to define the macro PARTICLE_SIMULATOR_TWO_DIMENSION)"
    parser.add_argument('-DPARTICLE_SIMULATOR_TWO_DIMENSION', \
                        required=False, help=help_msg, \
                        dest='PARTICLE_SIMULATOR_TWO_DIMENSION',\
                        action='store_true')
    # Return the result 
    args = parser.parse_args()
    if (args.PARTICLE_SIMULATOR_TWO_DIMENSION == False):
       dim_num = 3
    else:
       dim_num = 2
    return args.input_files, args.output_dir, dim_num

    # [DEBUG] Check
    #print("input is {}".format(args.input_files))
    #print("output_dir is {}".format(args.output_dir))
#================================
#   Main function
#================================
if __name__ == '__main__':
    # Analyze the command-line options and arguments
    input_files, output_dir, dim_num = analyze_CL_args()
    #print("input is {}".format(input_files))
    #print("output_dir is {}".format(output_dir))

    # Create auto-generator object
    gen = Automatic_Generator(input_files,output_dir,dim_num)
    gen.read_files(input_files)
    gen.analyze()
    gen.check()
    gen.generate(output_dir)
