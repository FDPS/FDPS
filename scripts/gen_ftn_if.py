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
    print("Module os,sys,struct,math are not found.")
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

#================================
#   Class definition
#================================
class Method_Generation_Data:
    def __init__(self):
        self.class_name   = ""
        self.member_pairs = []


class File_Data:
    def __init__(self):
        self.name    = ""
        self.content = []
        self.modules = []


class Module_Data:
    def __init__(self):
        self.name         = ""
        self.buf_addr_sta = -1
        self.buf_addr_end = -1
        self.structures   = []


class Structure_Data:
    def __init__(self):
        self.name         = ""
        self.attrib       = {"FP":False, \
                             "EPI":False, \
                             "EPJ":False, \
                             "Force":False}
        # Method database
        self.method_DB    = {"getId":[False,False,[]], \
                             "getPos":[False,False,[]], \
                             "setPos":[False,False,[]], \
                             "getCharge":[False,False,[]], \
                             "getChargePM":[False,False,[]], \
                             "getRSearch":[False,False,[]], \
                             "copyFromForce":[False,False,[]], \
                             "copyFromForcePM":[False,False,""], \
                             "copyFromFP":[False,False,[]], \
                             "clear":[False,False,[]]} 
        self._INDX_GENERABLE_FLAG = 0
        self._INDX_NECESSARY_FLAG = 1
        self._INDX_DATA_TO_GEN    = 2
        # where the 1st values are the values of the `generable` flags,
        # the 2nd values are the values of the `necessary` flags,
        # the 3rd values are the data required to generate the methods,
        # which correspond to the dummy arguments of private functions 
        # __write_* (where * is the name of method).
        # Note that each field of value should be accessed via 
        # self._INDX_*

        self.buf_addr_sta = -1
        self.buf_addr_end = -1
        self.members      = []

class Member_Data:
    def __init__(self):
        self.name       = ""
        self.attrib     = ""
        self.buf_addr   = -1
        self.data_type  = ""
        self.is_array   = False
        self.array_dim  = -1


class Automatic_Generator:
    def __init__(self,input_files,output_dir,dim_num):
        #=======================================================================
        # Display copyright
        text = """
        ################################################

        Python script to generate FDPS Fortran interface 
        (C) Copyright 2016  FDPS developer team. 

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
        self.__ftnDT_to_cppDT = {"integer(kind=c_int)":"int", \
                                 "integer(c_int)":"int", \
                                 "integer(kind=c_short)":"short int", \
                                 "integer(c_short)":"short int", \
                                 "integer(kind=c_long)":"long int", \
                                 "integer(c_long)":"long int", \
                                 "integer(kind=c_long_long)":"PS::S64", \
                                 "integer(c_long_long)":"PS::S64", \
                                 "integer(kind=c_signed_char)":"signed char", \
                                 "integer(c_signed_char)":"signed char", \
                                 "integer(kind=c_size_t)":"size_t", \
                                 "integer(c_size_t)":"size_t", \
                                 "integer(kind=c_int8_t)":"int8_t", \
                                 "integer(c_int8_t)":"int8_t", \
                                 "integer(kind=c_int16_t)":"int16_t", \
                                 "integer(c_int16_t)":"int16_t", \
                                 "integer(kind=c_int32_t)":"PS::S32", \
                                 "integer(c_int32_t)":"PS::S32", \
                                 "integer(kind=c_int64_t)":"PS::S64", \
                                 "integer(c_int64_t)":"PS::S64", \
                                 "integer(kind=c_int_least8_t)":"int_least8_t", \
                                 "integer(c_int_least8_t)":"int_least8_t", \
                                 "integer(kind=c_int_least16_t)":"int_least16_t", \
                                 "integer(c_int_least16_t)":"int_least16_t", \
                                 "integer(kind=c_int_least32_t)":"int_least32_t", \
                                 "integer(c_int_least32_t)":"int_least32_t", \
                                 "integer(kind=c_int_least64_t)":"int_least64_t", \
                                 "integer(c_int_least64_t)":"int_least64_t", \
                                 "integer(kind=c_int_fast8_t)":"int_fast8_t", \
                                 "integer(c_int_fast8_t)":"int_fast8_t", \
                                 "integer(kind=c_int_fast16_t)":"int_fast16_t", \
                                 "integer(c_int_fast16_t)":"int_fast16_t", \
                                 "integer(kind=c_int_fast32_t)":"int_fast32_t", \
                                 "integer(c_int_fast32_t)":"int_fast32_t", \
                                 "integer(kind=c_int_fast64_t)":"int_fast64_t", \
                                 "integer(c_int_fast64_t)":"int_fast64_t", \
                                 "integer(kind=c_intmax_t)":"intmax_t", \
                                 "integer(c_intmax_t)":"intmax_t", \
                                 "integer(kind=c_intptr_t)":"intptr_t", \
                                 "integer(c_intptr_t)":"intptr_t", \
                                 "real(kind=c_float)":"PS::F32", \
                                 "real(c_float)":"PS::F32", \
                                 "real(kind=c_double)":"PS::F64", \
                                 "real(c_double)":"PS::F64", \
                                 "real(kind=c_long_double)":"long double", \
                                 "real(c_long_double)":"long double", \
                                 "complex(kind=c_float_complex)":"float _Complex", \
                                 "complex(c_float_complex)":"float _Complex", \
                                 "complex(kind=c_double_complex)":"double _Complex", \
                                 "complex(c_double_complex)":"double _Complex", \
                                 "complex(kind=c_long_double_complex)":"long double _Complex", \
                                 "complex(c_long_double_complex)":"long double _Complex", \
                                 "logical(kind=c_bool)":"_Bool", \
                                 "logical(c_bool)":"_Bool", \
                                 "character(kind=c_char)":"char", \
                                 "character(c_char)":"char", \
                                 "type(fdps_f32vec)":"PS::F32vec", \
                                 "type(fdps_f64vec)":"PS::F64vec", \
                                 "type(fdps_f32mat)":"PS::F32mat", \
                                 "type(fdps_f64mat)":"PS::F64mat"}
        self.__usable_data_types = frozenset(self.__ftnDT_to_cppDT.keys())
        self.__default_IV_of_cppDT = {"int":0, \
                                      "short int":0, \
                                      "long int":0, \
                                      "long long int":0, \
                                      "PS::S32":0, \
                                      "PS::S64":0, \
                                      "size_t":0, \
                                      "int8_t":0, \
                                      "int16_t":0, \
                                      "int32_t":0, \
                                      "int64_t":0, \
                                      "int_least8_t":0, \
                                      "int_least16_t":0, \
                                      "int_least32_t":0, \
                                      "int_least64_t":0, \
                                      "int_fast8_t":0, \
                                      "int_fast16_t":0, \
                                      "int_fast32_t":0, \
                                      "int_fast64_t":0, \
                                      "intmax_t":0, \
                                      "intptr_t":0, \
                                      "float":0.0, \
                                      "double":0.0, \
                                      "long double":0.0, \
                                      "PS::F32":0.0, \
                                      "PS::F64":0.0, \
                                      "float _Complex":0.0, \
                                      "double _Complex":0.0, \
                                      "long double _Complex":0.0, \
                                      "_Bool":"false", \
                                      "char":"", \
                                      "PS::F32vec":0.0, \
                                      "PS::F64vec":0.0, \
                                      "PS::F32mat":0.0, \
                                      "PS::F64mat":0.0}
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
        # Fortran file manager
        self.files      = []
        #=======================================================================
        # Lists of FullParticle, EssentialParticleI, EssentialParticleJ, Force
        self.__FPs    = []
        self.__EPIs   = []
        self.__EPJs   = []
        self.__Forces = []
        # A list of SuperParticleJ classes
        self.__SPJs_ftn = ["fdps_spj_monopole", \
                           "fdps_spj_quadrupole", \
                           "fdps_spj_monopole_geomcen", \
                           "fdps_spj_dipole_geomcen", \
                           "fdps_spj_quadrupole_geomcen", \
                           "fdps_spj_monopole_scatter", \
                           "fdps_spj_quadrupole_scatter", \
                           "fdps_spj_monopole_cutoff"]
        self.__SPJs_cpp = ["PS::SPJMonopole", \
                           "PS::SPJQuadrupole", \
                           "PS::SPJMonopoleGeometricCenter", \
                           "PS::SPJDipoleGeometricCenter", \
                           "PS::SPJQuadrupoleGeometricCenter", \
                           "PS::SPJMonopoleScatter", \
                           "PS::SPJQuadrupoleScatter", \
                           "PS::SPJMonopoleCutoff"]
        # A list of the modes of tree with PS::SEARCH_MODE_LONG*
        self.__multipole_kinds = ["Monopole", \
                                  "Quadrupole", \
                                  "MonopoleGeometricCenter", \
                                  "DipoleGeometricCenter", \
                                  "QuadrupoleGeometricCenter", \
                                  "MonopoleWithScatterSearch", \
                                  "QuadrupoleWithScatterSearch", \
                                  "MonopoleWithCutoff"]
        # A list of the modes of tree with short-range interaction
        self.__neighbor_kinds = ["Gather", \
                                 "Scatter", \
                                 "Symmetry"]
        # A list of tree types
        self.__tree_kinds = []
        # Lists of calcForce*()
        self.__calc_force_kinds_short = []
        self.__calc_force_kinds_long  = []

        # List of data types supported by reduction or comminucation functions
        self.__comm_int_types = collections.OrderedDict()
        self.__comm_int_types["integer(kind=c_int)"]        = "i32"
        self.__comm_int_types["integer(kind=c_long_long)"]  = "i64"
        self.__comm_real_types = collections.OrderedDict()
        self.__comm_real_types["real(kind=c_float)"]        = "r32"
        self.__comm_real_types["real(kind=c_double)"]       = "r64"
        self.__comm_all_types = collections.OrderedDict()
        self.__comm_all_types.update(self.__comm_int_types)
        self.__comm_all_types.update(self.__comm_real_types)

    def __print_warning(self,msg):
        #output_msg = "\033[33;5m" + "[warning] " + msg + "\033[m" # yellow blink
        output_msg = "\033[33;1m" + "[warning] " + msg + "\033[m" # yellow
        print(output_msg)

    def __print_error(self,msg):
        #output_msg = "\033[31;5m" + "[error] " + msg + "\033[m" # red blink
        output_msg = "\033[31;1m" + "[error] " + msg + "\033[m" # red
        print(output_msg)

    def __print_checkpoint(self,msg):
        #output_msg = "\033[32;5m" + "[check point] " + msg + "\033[m" # green blink
        output_msg = "\033[32;1m" + "[check point] " + msg + "\033[m" # green
        print(output_msg)

    def __print_dbgmsg(self,msg):
        #output_msg = "\033[35;5m" + "[debug] " + msg + "\033[m" # magenta blink
        output_msg = "\033[35;1m" + "[debug] " + msg + "\033[m" # magenta
        print(output_msg)


    def __module_existence_check(self,module_name):
        bool_val = False
        for f in self.files:
            for mod in f.modules:
                if (module_name == mod.name):
                    bool_val = True
        return bool_val

    def __class_existence_check(self,class_name):
        bool_val = False
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (class_name == s.name):
                        bool_val = True
        return bool_val 

    def __member_existence_check(self,class_name,member_name):
        bool_val = False
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (class_name == s.name):
                        for mbr in s.members:
                            if (member_name == mbr.name):
                                bool_val = True
        return bool_val

    def __getId_existence_check(self,class_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (class_name == s.name):
                        # First, we check if class_name is EPJ or not
                        if (s.attrib["EPJ"] == False):
                            msg = "{0} is not EPJ!".format(class_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            if (s.method_DB["getId"][s._INDX_GENERABLE_FLAG]):
                                return True
                            else:
                                return False

    def __getRSearch_existence_check(self,class_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (class_name == s.name):
                        # First, we check if class_name is either EPJ or FP
                        if ((s.attrib["FP"] == False) and \
                            (s.attrib["EPJ"] == False)):
                            msg = "{0} is neither EPJ nor FP!".format(class_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            if (s.method_DB["getRSearch"][s._INDX_GENERABLE_FLAG]):
                                return True
                            else:
                                return False

    def __check_copyFromFP_consistency(self,EP_name,FP_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (EP_name == s.name):
                        if ((s.attrib["EPI"] == False) and \
                            (s.attrib["EPJ"] == False)):
                            msg = "{0} is neither EPI nor EPJ!".format(EP_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            if (s.method_DB["copyFromFP"][s._INDX_GENERABLE_FLAG]):
                                meth_gen_data = s.method_DB["copyFromFP"][s._INDX_DATA_TO_GEN]
                                ret = False
                                for data in meth_gen_data:
                                    class_name = data.class_name
                                    if (FP_name == class_name):
                                         ret = True
                                return ret
                            else:
                                return False

    def __check_copyFromForce_consistency(self,FP_name,Force_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (FP_name == s.name):
                        if (s.attrib["FP"] == False):
                            msg = "{0} is not FP!".format(FP_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            if (s.method_DB["copyFromForce"][s._INDX_GENERABLE_FLAG]):
                                meth_gen_data = s.method_DB["copyFromForce"][s._INDX_DATA_TO_GEN]
                                ret = False
                                for data in meth_gen_data:
                                    class_name = data.class_name
                                    if (Force_name == class_name):
                                        ret = True
                                return ret
                            else:
                                return False

    def __get_mod_name_by_FP(self,FP_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (FP_name == s.name):
                        if (s.attrib["FP"] == False):
                            msg = "{0} is not FP!".format(FP_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            return mod.name

    def __get_mod_name_by_EPI(self,EPI_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (EPI_name == s.name):
                        if (s.attrib["EPI"] == False):
                            msg = "{0} is not EPI!".format(EPI_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            return mod.name

    def __get_mod_name_by_EPJ(self,EPJ_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (EPJ_name == s.name):
                        if (s.attrib["EPJ"] == False):
                            msg = "{0} is not EPJ!".format(EPJ_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            return mod.name

    def __get_mod_name_by_Force(self,Force_name):
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    if (Force_name == s.name):
                        if (s.attrib["Force"] == False):
                            msg = "{0} is not Force!".format(Force_name)
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            return mod.name

    def __is_s32(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "integer(kind=c_int)") and \
                          (mbr.is_array == False))
        conditions.append((mbr.data_type == "integer(c_int)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_s64(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "integer(kind=c_long_long)") and \
                          (mbr.is_array == False))
        conditions.append((mbr.data_type == "integer(c_long_long)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f32(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "real(kind=c_float)") and \
                          (mbr.is_array == False))
        conditions.append((mbr.data_type == "real(c_float)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f64(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "real(kind=c_double)") and \
                          (mbr.is_array == False))
        conditions.append((mbr.data_type == "real(c_double)") and \
                          (mbr.is_array == False))
        return any(conditions)

    def __is_f32vec(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "type(fdps_f32vec)"))
        conditions.append((mbr.data_type == "real(kind=c_float)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        conditions.append((mbr.data_type == "real(c_float)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        return any(conditions)

    def __is_f64vec(self,mbr):
        conditions = []
        conditions.append((mbr.data_type == "type(fdps_f64vec)"))
        conditions.append((mbr.data_type == "real(kind=c_double)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        conditions.append((mbr.data_type == "real(c_double)") and \
                           mbr.is_array and \
                          (mbr.array_dim == self.__DIM_NUM))
        return any(conditions)

    def __get_fdpsDT(self,mbr):
        # This method converts Fortran data type to data type used in FDPS.
        # This method works correctly only if data in `mbr` is complete.
        if (mbr.is_array):
            if (mbr.attrib not in self.__usable_fdps_mbr_dirs):
                # normal array cases
                data_type = self.__ftnDT_to_cppDT[mbr.data_type]
                array_dim = mbr.array_dim
            else:
                # `position` or `velocity`
                if (self.__is_f32vec(mbr)):
                    data_type = "PS::F32vec"
                    array_dim = 0
                else:
                    data_type = "PS::F64vec"
                    array_dim = 0
        else:
            data_type = self.__ftnDT_to_cppDT[mbr.data_type]
            array_dim = 0
        return data_type,array_dim

    def __get_ftn_indent(self,level):
        text = "   " * int(level)
        return text

    def __get_cpp_indent(self,level):
        text = "   " * int(level)
        return text

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

    def __write_getId(self,fh,data_type,member_name):
        text = """
        {DATA_TYPE} getId() const {{ 
           return this->{MBR_NAME};
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)

    def __write_getPos(self,fh,data_type,member_name):
        text = """
        {DATA_TYPE} getPos() const {{ 
           return this->{MBR_NAME};
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)
    
    def __write_setPos(self,fh,data_type,member_name):
        text = """
        void setPos(const {DATA_TYPE} pos_new) {{
           this->{MBR_NAME} = pos_new;
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)
    
    def __write_getCharge(self,fh,data_type,member_name):
        text = """
        {DATA_TYPE} getCharge() const {{
           return this->{MBR_NAME};
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)
    
    def __write_getChargePM(self,fh,data_type,member_name):
        text = """
        {DATA_TYPE} getChargeParticleMesh() const {{
           return this->{MBR_NAME};
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)
       
    def __write_getRSearch(self,fh,data_type,member_name):
        text = """
        {DATA_TYPE} getRSearch() const {{
           return this->{MBR_NAME};
        }}
        """.format(DATA_TYPE=data_type, \
                   MBR_NAME=member_name)
        text = text[1:].rstrip()
        text = textwrap.dedent(text)
        indent = self.__get_cpp_indent(1)
        text = self.__add_indent(text,indent)
        fh.write(text)
    
    def __write_copyFromForce(self,fh,meth_gen_data):
        text = ""
        for data in meth_gen_data:
           class_name   = data.class_name
           member_pairs = data.member_pairs
           text += "   void copyFromForce(const {0} & force) {{\n".format(class_name)
           for items in member_pairs:
               text += "      this->{0} = force.{1};\n".format(items[0],items[1])
           text += "   }\n"
        fh.write(text)
    
    def __write_copyFromForcePM(self,fh,member_name):
        text = "   void copyFromForceParticleMesh(const PS::F32vec & acc_pm) {\n" \
             + "      this->{0} = acc_pm;\n".format(member_name) \
             + "   }\n"
        fh.write(text)
    
    def __write_copyFromFP(self,fh,meth_gen_data):
        text = ""
        for data in meth_gen_data:
            class_name   = data.class_name
            member_pairs = data.member_pairs
            text += "   void copyFromFP(const {0} & fp) {{\n".format(class_name) 
            for items in member_pairs:
                text += "      this->{0} = fp.{1};\n".format(items[0],items[1]) 
            text += "   }\n"
        fh.write(text)
    
    def __write_clear(self,fh,meth_gen_data):
        # Extract gen_mode
        gen_mode = meth_gen_data[0]
        # Branches by gen_mode
        if (gen_mode == 0):
            # In these cases, we generate a normal function.
            # (i) Get data to generate
            members    = meth_gen_data[1]
            data_types = meth_gen_data[2]
            array_dims = meth_gen_data[3]
            vals       = meth_gen_data[4]
            num_mbr    = len(members)
            # (ii) Set text
            text = "   void clear() {\n"
            for i in range(0,num_mbr):
                mbr_name  = members[i]
                data_type = data_types[i]
                array_dim = array_dims[i]
                val       = vals[i]
                if (array_dim == 0):
                    # Scalar case
                    text += "      this->{0} = {1};\n".format(mbr_name,val)
                else:
                    # Array case
                    if (data_type != "char"):
                        text += "      for (int i=0; i<{0}; i++)\n".format(array_dim)
                        text += "         this->{0}[i] = {1};\n".format(mbr_name,val)
                    else:
                        text += "      std::strcpy({0},{1});\n".format(mbr_name,val)
            text += "   }\n"
        elif (gen_mode == 1):
            # In this case, we generate the declaration of clear method.
            # The implementation of the clear method generated by
            # a different function in the Python script.
            text = "   void clear();\n"
        fh.write(text)

    def __write_clear_impl(self,fh,meth_gen_data,force_name):
        # Extract gen_mode
        gen_mode = meth_gen_data[0]
        # Generate the implementation of the clear method.
        if (gen_mode == 1):
            subroutine_name = meth_gen_data[1]
            text = """
            // Impl. of the clear method of the class `{FORCE}`
            extern "C" void {FTN_SUBROUTINE}({FORCE} *f);
            void {FORCE}::clear() {{
                {FTN_SUBROUTINE}(this);
            }};
            """.format(FTN_SUBROUTINE=subroutine_name, \
                       FORCE=force_name)
            text = text[1:].rstrip() + "\n\n"
            text = textwrap.dedent(text)
            fh.write(text)

    def __write_psys_branch_for_create_psys(self,fh):
        # Set the indents
        indent = self.__get_cpp_indent(2)
        # Set the output text
        texts = []
        is_first = True
        for fp_t in self.__FPs:
            if (is_first):
                if_chars = "if"
                is_first = False
            else:
                if_chars = "} else if"
            # if & elif clauses
            text = """
            {IF_CHARS} (psys_info_ == "{FP_T}") {{
               PS::ParticleSystem<{FP_T}> *psys;
               psys = new PS::ParticleSystem<{FP_T}>;
               psys_data.id   = num_psys_creation;
               psys_data.ptr  = (void *) psys;
               psys_data.info = psys_info_;
            """.format(IF_CHARS=if_chars, \
                       FP_T=fp_t)
            text = text[1:].rstrip() + "\n" 
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # else-clause
        text = """
        } else { 
           PS::S32 myrank = PS::Comm::getRank();
           if (myrank == 0) {
               std::string errmsg,func_name;
               errmsg = \"FullParticle `\" + psys_info_ + \"` does not exist.\";
               func_name = "create_psys";
               print_errmsg(errmsg,func_name);
           }
           PS::Abort(-1);
           std::exit(EXIT_FAILURE);
        }
        """[1:].rstrip() + "\n"
        text = textwrap.dedent(text)
        text = self.__add_indent(text,indent)
        texts.append(text)
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_psys_branch_for_normal_APIs(self,fh,inst,cmnt=""):
        # Set the indent
        indent = self.__get_cpp_indent(4)
        # Set the output text
        texts = []
        is_first = True
        for fp_t in self.__FPs:
            if (is_first):
                if_chars = "if"
                is_first = False
            else:
                if_chars = "} else if"
            # (1) if & elif clauses
            text = """
            {IF_CHARS} (it->info == "{FP_T}") {{
               PS::ParticleSystem<{FP_T}> *psys;
               psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
               {INSTRUCTION}
               {COMMENT}
            """.format(IF_CHARS=if_chars, \
                       FP_T=fp_t, \
                       INSTRUCTION=inst, \
                       COMMENT=cmnt)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) else clause
        text = """
        } else {
           PS::S32 myrank = PS::Comm::getRank();
           if (myrank == 0) {
              std::string errmsg,func_name;
              errmsg = "An invalid ParticleSystem number is received.";
              func_name = "";
              print_errmsg(errmsg,func_name);
           }
           PS::Finalize();
           std::exit(EXIT_FAILURE);
        }
        break;
        """[1:].rstrip() + "\n"
        text = textwrap.dedent(text)
        text = self.__add_indent(text,indent)
        texts.append(text)
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_cpp_impl(self,fh):
        # Set the common indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "add_particle"
        func_num = 0
        for fp_t in self.__FPs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int psys_num,
            {INDENT}const {FP_T} *ptcl) {{
               for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                  if (it->id == psys_num) {{
                     PS::ParticleSystem<{FP_T}> *psys;
                     psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                     psys->addOneParticle(*ptcl);
                     break;
                  }}
               }}
            }}
            """.format(FUNC_NAME=func_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_sort_particle_cpp_impl(self,fh):
        # Set the common indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "sort_particle"
        func_num = 0
        for fp_t in self.__FPs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int psys_num,
            {INDENT}bool (*pfunc_comp)(const struct {FP_T} &,
            {INDENT}                   const struct {FP_T} &)) {{
               for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                  if (it->id == psys_num) {{
                     PS::ParticleSystem<{FP_T}> *psys;
                     psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                     psys->sortParticle(pfunc_comp);
                     break;
                  }}
               }}
            }}
            """.format(FUNC_NAME=func_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_tree_branch_for_create_tree(self,fh):
        # Set the indent
        indent = self.__get_cpp_indent(2)
        # Set the output text
        texts = []
        is_first = True
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            if (is_first):
                if_chars = "if"
                is_first = False
            else:
                if_chars = "} else if"
            # (1) if & elif clauses
            text = """
            {IF_CHARS} (tree_info_ == "{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MODE_T}") {{
               PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *tree;
               tree = new PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T};
               tree_data.id   = num_tree_creation;
               tree_data.ptr  = (void *) tree;
               tree_data.info = tree_info_;
            """.format(IF_CHARS=if_chars, \
                       TREE_T=tree_t, \
                       FORCE_T=force_t, \
                       EPI_T=epi_t, \
                       EPJ_T=epj_t, \
                       MODE_T=mode_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) else clause
        text = """
        } else {
           PS::S32 myrank = PS::Comm::getRank();
           if (myrank == 0) {
              std::string errmsg,func_name;
              errmsg = "cannot create Tree `" + tree_info_ + "` ";
              func_name = "create_tree";
              print_errmsg(errmsg,func_name);
           }
           PS::Finalize();
           std::exit(EXIT_FAILURE);
        }
        """[1:].rstrip() + "\n"
        text = textwrap.dedent(text)
        text = self.__add_indent(text,indent)
        texts.append(text)
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_tree_branch_for_normal_APIs(self,fh,inst):
        # Set the indent
        indent = self.__get_cpp_indent(4)
        # Set the output text
        texts = []
        is_first = True
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            if (is_first):
                if_chars = "if"
                is_first = False
            else:
                if_chars = "} else if"
            # (1) if & elif clauses
            text = """
            {IF_CHARS} (it->info == "{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MODE_T}") {{
               PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *tree;
               tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *) it->ptr;
               {INSTRUCTION}
            """.format(IF_CHARS=if_chars, \
                       TREE_T=tree_t, \
                       FORCE_T=force_t, \
                       EPI_T=epi_t, \
                       EPJ_T=epj_t, \
                       MODE_T=mode_t, \
                       INSTRUCTION=inst)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) else clause
        text = """
        } else {
           PS::S32 myrank = PS::Comm::getRank();
           if (myrank == 0) {
              std::string errmsg,func_name;
              errmsg = "An invalid Tree number is received.";
              func_name = "";
              print_errmsg(errmsg,func_name);
           }
           PS::Abort(-1);
           std::exit(EXIT_FAILURE);
        }
        """[1:].rstrip() + "\n"
        text = textwrap.dedent(text)
        text = self.__add_indent(text,indent)
        texts.append(text)
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_branch(self,fh,generic_func_name,fdps_method_name):
        # Check arguments
        ret = isinstance(generic_func_name,str) and \
              isinstance(fdps_method_name,str)
        if (ret == False):
           print("[dev] invalid arguments in __write_calc_force_branch()!")
           sys.exit()
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        #=================================
        # (1) Short-range interaction
        #=================================
        func_name_base = generic_func_name + "_s"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear,
                {INDENT}const enum PS_INTERACTION_LIST_MODE list_mode) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get psys
                   PS::ParticleSystem<{FP_T}> *psys;
                   for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                      if (it->id == psys_num) {{
                         psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                         break;
                      }}
                   }}
                   // Get dinfo
                   PS::DomainInfo *dinfo;
                   for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {{
                      if (it->id == dinfo_num) {{
                         dinfo = (PS::DomainInfo *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   if (list_mode == MAKE_LIST) {{
                      tree->{API_NAME}(pfunc_ep_ep, *psys, *dinfo, clear, PS::MAKE_LIST);
                   }} else if (list_mode == MAKE_LIST_FOR_REUSE) {{
                      tree->{API_NAME}(pfunc_ep_ep, *psys, *dinfo, clear, PS::MAKE_LIST_FOR_REUSE);
                   }} else if (list_mode == REUSE_LIST) {{
                      tree->{API_NAME}(pfunc_ep_ep, *psys, *dinfo, clear, PS::REUSE_LIST);
                   }} else {{
                      PS::S32 myrank = PS::Comm::getRank();
                      if (myrank == 0) {{
                         std::string errmsg,func_name;
                         errmsg = "Unknow list_mode is specified.";
                         func_name = "{FUNC_NAME}";
                         print_errmsg(errmsg,func_name);
                      }}
                      PS::Abort(-1);
                      std::exit(EXIT_FAILURE);
                   }}
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1;
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const bool clear) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get psys
                   PS::ParticleSystem<{FP_T}> *psys;
                   for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                      if (it->id == psys_num) {{
                         psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   tree->{API_NAME}(pfunc_ep_ep, *psys, clear);
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                # update func_num
                func_num += 1;
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,search_t in self.__tree_kinds:
                if (tree_t == "Long"):
                    continue
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{SEARCH_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get dinfo
                   PS::DomainInfo *dinfo;
                   for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {{
                      if (it->id == dinfo_num) {{
                         dinfo = (PS::DomainInfo *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   tree->{API_NAME}(pfunc_ep_ep, *dinfo, clear);
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1;

        #=================================
        # (2) Long-range interaction
        #=================================
        func_name_base = generic_func_name + "_l"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear,
                {INDENT}const enum PS_INTERACTION_LIST_MODE list_mode) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get psys
                   PS::ParticleSystem<{FP_T}> *psys;
                   for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                      if (it->id == psys_num) {{
                         psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                         break;
                      }}
                   }}
                   // Get dinfo
                   PS::DomainInfo *dinfo;
                   for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {{
                      if (it->id == dinfo_num) {{
                         dinfo = (PS::DomainInfo *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   if (list_mode == MAKE_LIST) {{
                      tree->{API_NAME}(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo, clear, PS::MAKE_LIST);
                   }} else if (list_mode == MAKE_LIST_FOR_REUSE) {{
                      tree->{API_NAME}(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo, clear, PS::MAKE_LIST_FOR_REUSE);
                   }} else if (list_mode == REUSE_LIST) {{
                      tree->{API_NAME}(pfunc_ep_ep, pfunc_ep_sp, *psys, *dinfo, clear, PS::REUSE_LIST);
                   }} else {{
                      PS::S32 myrank = PS::Comm::getRank();
                      if (myrank == 0) {{
                         std::string errmsg,func_name;
                         errmsg = "Unknow list_mode is specified.";
                         func_name = "{FUNC_NAME}";
                         print_errmsg(errmsg,func_name);
                      }}
                      PS::Abort(-1);
                      std::exit(EXIT_FAILURE);
                   }}
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const bool clear) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get psys
                   PS::ParticleSystem<{FP_T}> *psys;
                   for (std::vector<Psys_Data>::iterator it = psys_vector.begin(); it != psys_vector.end(); ++it) {{
                      if (it->id == psys_num) {{
                         psys = (PS::ParticleSystem<{FP_T}> *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   tree->{API_NAME}(pfunc_ep_ep, pfunc_ep_sp, *psys, clear);
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,multipole_t in self.__tree_kinds:
                if (tree_t == "Short"):
                    continue
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear) {{
                   // Get tree
                   PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *tree;
                   for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                      if (it->id == tree_num) {{
                         tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MULTIPOLE_T} *) it->ptr;
                         break;
                      }}
                   }}
                   // Get dinfo
                   PS::DomainInfo *dinfo;
                   for (std::vector<Dinfo_Data>::iterator it = dinfo_vector.begin(); it != dinfo_vector.end(); ++it) {{
                      if (it->id == dinfo_num) {{
                         dinfo = (PS::DomainInfo *) it->ptr;
                         break;
                      }}
                   }}
                   // Issue calcForce*()
                   tree->{API_NAME}(pfunc_ep_ep, pfunc_ep_sp, *dinfo, clear);
                }}
                """.format(FUNC_NAME=func_name,\
                           API_NAME=fdps_method_name, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_branch(self,fh):
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "get_neighbor_list"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            if (tree_t == "Long"):
                if ((mode_t != "MonopoleWithScatterSearch") and \
                    (mode_t != "QuadrupoleWithScatterSearch")):
                    continue
            # Generate
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int tree_num,
            {INDENT}const PS::F64vec *pos,
            {INDENT}const PS::F64 r_search,
            {INDENT}int *num_epj,
            {INDENT}void **cptr_to_epj) {{
               // Get tree
               PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *tree;
               for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                  if (it->id == tree_num) {{
                     tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *) it->ptr;
                     break;
                  }}
               }}
               // Get neighbor list
               struct NBS_Target_Data ptcl;
               struct {EPJ_T} * epj;
               ptcl.pos = *pos;
               ptcl.r_search = r_search;
               *num_epj = tree->getNeighborListOneParticle(ptcl,epj);
               *cptr_to_epj = (void *) &(epj[0]);
            }}
            """.format(FUNC_NAME=func_name,\
                       TREE_T=tree_t, \
                       FORCE_T=force_t, \
                       EPI_T=epi_t, \
                       EPJ_T=epj_t, \
                       MODE_T=mode_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_branch(self,fh):
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "get_epj_from_id"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            # [Note] 
            #    If EPJ does not have getId(), we do not generate this API.
            if (not self.__getId_existence_check(epj_t)):
                continue
            # Generate
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int tree_num,
            {INDENT}const PS::S64 id,
            {INDENT}void **cptr_to_epj) {{
               // Get tree
               PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *tree;
               for (std::vector<Tree_Data>::iterator it = tree_vector.begin(); it != tree_vector.end(); ++it) {{
                  if (it->id == tree_num) {{
                     tree = (PS::TreeForForce{TREE_T}<{FORCE_T},{EPI_T},{EPJ_T}>::{MODE_T} *) it->ptr;
                     break;
                  }}
               }}
               // Get EPJ
               struct {EPJ_T} * epj;
               epj = tree->getEpjFromId(id);
               *cptr_to_epj = (void *) epj;
            }}
            """.format(FUNC_NAME=func_name,\
                       TREE_T=tree_t, \
                       FORCE_T=force_t, \
                       EPI_T=epi_t, \
                       EPJ_T=epj_t, \
                       MODE_T=mode_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_protoDCL(self,fh):
        # Set the common indent
        indent = self.__get_cpp_indent(1)
        # Set the output texts
        texts = []
        func_name_base = "add_particle"
        func_num = 0
        for fp_t in self.__FPs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int psys_num,
            {INDENT}const {FP_T} *ptcl);
            """.format(FUNC_NAME=func_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_sort_particle_protoDCL(self,fh):
        # Set the common indent
        indent = self.__get_cpp_indent(1)
        # Set the output texts
        texts = []
        func_name_base = "sort_particle"
        func_num = 0
        for fp_t in self.__FPs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6)
            text = """
            void {FUNC_NAME}(const int psys_num,
            {INDENT}bool (*pfunc_comp)(const struct {FP_T} & left,
            {INDENT}                   const struct {FP_T} & right));
            """.format(FUNC_NAME=func_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_protoDCL_branch(self,fh,generic_func_name):
        # Check arguments
        ret = isinstance(generic_func_name,str)
        if (ret == False):
           print("[dev] invalid arguments in __write_calc_force_protoDCL_branch()!")
           sys.exit()
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        #=================================
        # (1) Short-range interaction
        #=================================
        func_name_base = generic_func_name + "_s"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear,
                {INDENT}const enum PS_INTERACTION_LIST_MODE list_mode);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1;
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const bool clear);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1;
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,search_t in self.__tree_kinds:
                if (tree_t == "Long"):
                    continue
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1;


        #=================================
        # (2) Long-range interaction
        #=================================
        func_name_base = generic_func_name + "_l"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear,
                {INDENT}const enum PS_INTERACTION_LIST_MODE list_mode);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int psys_num,
                {INDENT}const bool clear);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,multipole_t in self.__tree_kinds:
                if (tree_t == "Short"):
                    continue
                func_name = func_name_base + "{0:05d}".format(func_num)
                indent_at_func_args = " " * int(len(func_name) + 6 + 7)
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                text = """
                extern void {FUNC_NAME}(const int tree_num,
                {INDENT}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {EPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT}                    int ,
                {INDENT}                    {SPJ_T} *,
                {INDENT}                    int ,
                {INDENT}                    struct {FORCE_T} *),
                {INDENT}const int dinfo_num,
                {INDENT}const bool clear);
                """.format(FUNC_NAME=func_name,\
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_protoDCL_branch(self,fh):
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "get_neighbor_list"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            if (tree_t == "Long"):
                if ((mode_t != "MonopoleWithScatterSearch") and \
                    (mode_t != "QuadrupoleWithScatterSearch")):
                    continue
            # Generate
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6 + 7)
            text = """
            extern void {FUNC_NAME}(const int tree_num,
            {INDENT}const PS::F64vec *pos,
            {INDENT}const PS::F64 r_search,
            {INDENT}int *num_epj,
            {INDENT}void **cptr_to_epj);
            """.format(FUNC_NAME=func_name,\
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_protoDCL_branch(self,fh):
        # Set the indent
        indent = self.__get_cpp_indent(1)
        # Set the output text
        texts = []
        func_name_base = "get_epj_from_id"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            # [Note] 
            #    If EPJ does not have getId(), we do not generate this API.
            if (not self.__getId_existence_check(epj_t)):
                continue
            # Generate
            func_name = func_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(func_name) + 6 + 7)
            text = """
            extern void {FUNC_NAME}(const int tree_num,
            {INDENT}const PS::S64 id,
            {INDENT}void **cptr_to_epj);
            """.format(FUNC_NAME=func_name,\
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_cpp_if(self,fh):
        # Set the common indent
        indent = ""
        # Set the output texts
        texts = []
        real_name_base = "add_particle"
        func_num = 0
        for fp_t in self.__FPs:
            real_name = real_name_base + "{0:05d}".format(func_num)
            if_name  = "fdps_" + real_name
            indent_at_func_args = " " * int(len(if_name) + 6)
            text = """
            void {IFC_NAME}(const int psys_num,
            {INDENT}const {FP_T} *ptcl) {{
               {REAL_NAME}(psys_num,ptcl);
            }}
            """.format(IFC_NAME=if_name, \
                       REAL_NAME=real_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_sort_particle_cpp_if(self,fh):
        # Set the common indent
        indent = ""
        # Set the output texts
        texts = []
        real_name_base = "sort_particle"
        func_num = 0
        for fp_t in self.__FPs:
            real_name = real_name_base + "{0:05d}".format(func_num)
            if_name  = "fdps_" + real_name
            indent_at_func_args = " " * int(len(if_name) + 6)
            text = """
            void {IFC_NAME}(const int psys_num,
            {INDENT}bool (*pfunc_comp)(const struct {FP_T} & left,
            {INDENT}                   const struct {FP_T} & right)) {{
               {REAL_NAME}(psys_num,pfunc_comp);
            }}
            """.format(IFC_NAME=if_name, \
                       REAL_NAME=real_name, \
                       FP_T=fp_t, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_cpp_if(self,fh,generic_func_name):
        # Check arguments
        ret = isinstance(generic_func_name,str)
        if (ret == False):
           print("[dev] invalid arguments in __write_calc_force_cpp_if()!")
           sys.exit()
        # Set the indent
        indent = ""
        # Set the output texts
        texts = []
        #=================================
        # (1) Short-range interaction
        #=================================
        if_name_base = "fdps_" + generic_func_name + "_s"
        real_name_base = generic_func_name + "_s"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set the output text
                text = """
                void {IFC_NAME}(const int tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int psys_num,
                {INDENT1}const int dinfo_num,
                {INDENT1}const bool clear,
                {INDENT1}const enum PS_INTERACTION_LIST_MODE list_mode) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}psys_num,
                {INDENT2}dinfo_num,
                {INDENT2}clear,
                {INDENT2}list_mode);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set the output text
                text = """
                void {IFC_NAME}(const int tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int psys_num,
                {INDENT1}const bool clear) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}psys_num,
                {INDENT2}clear);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,search_t in self.__tree_kinds:
                if (tree_t == "Long"):
                    continue
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set the output text
                text = """
                void {IFC_NAME}(const int tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int dinfo_num,
                {INDENT1}const bool clear) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}dinfo_num,
                {INDENT2}clear);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        #=================================
        # (2) Long-range interaction
        #=================================
        if_name_base = "fdps_" + generic_func_name + "_l"
        real_name_base = generic_func_name + "_l"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set SuperParticleJ type
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                # Set the output text
                text = """
                void {IFC_NAME}(const int tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    {SPJ_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int psys_num,
                {INDENT1}const int dinfo_num,
                {INDENT1}const bool clear,
                {INDENT1}const enum PS_INTERACTION_LIST_MODE list_mode) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}pfunc_ep_sp,
                {INDENT2}psys_num,
                {INDENT2}dinfo_num,
                {INDENT2}clear,
                {INDENT2}list_mode);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set SuperParticleJ type
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                # Set the output text
                text = """
                void {IFC_NAME}(const int tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    {SPJ_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int psys_num,
                {INDENT1}const bool clear) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}pfunc_ep_sp,
                {INDENT2}psys_num,
                {INDENT2}clear);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,multipole_t in self.__tree_kinds:
                if (tree_t == "Short"):
                    continue
                # Set function names
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                real_name = real_name_base + func_label
                # Set the other indents
                indent_at_func_args1 = " " * int(len(if_name) + 6)
                indent_at_func_args2 = " " * int(len(real_name) + 4)
                # Set SuperParticleJ type
                spj_t = self.__SPJs_cpp[self.__multipole_kinds.index(multipole_t)]
                # Set the output text
                text = """
                void {IFC_NAME}(int const tree_num,
                {INDENT1}void (*pfunc_ep_ep)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {EPJ_T} *,
                {INDENT1}                    int , 
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}void (*pfunc_ep_sp)(struct {EPI_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    {SPJ_T} *,
                {INDENT1}                    int ,
                {INDENT1}                    struct {FORCE_T} *),
                {INDENT1}const int dinfo_num,
                {INDENT1}const bool clear) {{
                   {REAL_NAME}(tree_num, 
                {INDENT2}pfunc_ep_ep,
                {INDENT2}pfunc_ep_sp,
                {INDENT2}dinfo_num,
                {INDENT2}clear);
                }}
                """.format(IFC_NAME=if_name, \
                           REAL_NAME=real_name, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SPJ_T=spj_t, \
                           FORCE_T=force_t, \
                           INDENT1=indent_at_func_args1, \
                           INDENT2=indent_at_func_args2)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_cpp_if(self,fh):
        # Set the common indent
        indent = ""
        # Set the output texts
        texts = []
        generic_func_name = "get_neighbor_list"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            if (tree_t == "Long"):
                if ((mode_t != "MonopoleWithScatterSearch") and \
                    (mode_t != "QuadrupoleWithScatterSearch")):
                    continue
            # Generate
            real_name = generic_func_name + "{0:05d}".format(func_num)
            if_name  = "fdps_" + real_name
            indent_at_func_args1 = " " * int(len(if_name) + 6)
            indent_at_func_args2 = " " * int(len(real_name) + 4)
            text = """
            void {IFC_NAME}(const int tree_num,
            {INDENT1}const PS::F64vec *pos,
            {INDENT1}const PS::F64 r_search,
            {INDENT1}int *num_epj,
            {INDENT1}void **cptr_to_epj) {{
               {REAL_NAME}(tree_num,
            {INDENT2}pos,
            {INDENT2}r_search,
            {INDENT2}num_epj,
            {INDENT2}cptr_to_epj);
            }}
            """.format(IFC_NAME=if_name,\
                       REAL_NAME=real_name, \
                       INDENT1=indent_at_func_args1, \
                       INDENT2=indent_at_func_args2)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_cpp_if(self,fh):
        # Set the common indent
        indent = ""
        # Set the output texts
        texts = []
        generic_func_name = "get_epj_from_id"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            # [Note] 
            #    If EPJ does not have getId(), we do not generate this API.
            if (not self.__getId_existence_check(epj_t)):
                continue
            # Generate
            real_name = generic_func_name + "{0:05d}".format(func_num)
            if_name  = "fdps_" + real_name
            indent_at_func_args1 = " " * int(len(if_name) + 6)
            indent_at_func_args2 = " " * int(len(real_name) + 4)
            text = """
            void {IFC_NAME}(const int tree_num,
            {INDENT1}const PS::S64 id,
            {INDENT1}void **cptr_to_epj) {{
               {REAL_NAME}(tree_num,
            {INDENT2}id,
            {INDENT2}cptr_to_epj);
            }}
            """.format(IFC_NAME=if_name,\
                       REAL_NAME=real_name, \
                       INDENT1=indent_at_func_args1, \
                       INDENT2=indent_at_func_args2)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1;
        # Output
        text = "".join(texts)
        fh.write(text)
            

    def __write_reduction_cpp_if(self,fh, \
                                 generic_func_name, \
                                 fdps_api_name, \
                                 gen_func_w_id):
        # Set the common indent
        indent = ""
        # Set the output texts 
        texts = []
        # (1) Case w/o ID
        for key,label in self.__comm_all_types.items():
            data_type = self.__ftnDT_to_cppDT[key]
            func_name = "fdps_" + generic_func_name + "_" + label
            text = """
            void {FUNC_NAME}(const {DATA_TYPE} f_in, {DATA_TYPE} *f_out) {{
                //*f_out = {API_NAME}(f_in);
                {DATA_TYPE} tmp = f_in;
                *f_out = {API_NAME}(tmp);
            }}
            """.format(FUNC_NAME=func_name, \
                       DATA_TYPE=data_type, \
                       API_NAME=fdps_api_name)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) Case w/ ID
        if (gen_func_w_id):
            for key,label in self.__comm_real_types.items():
                data_type = self.__ftnDT_to_cppDT[key]
                func_name = "fdps_" + generic_func_name + "_w_id_" + label
                indent_at_args = " " * int(len(func_name) + 6)
                text = """
                void {FUNC_NAME}(const {DATA_TYPE} f_in,
                {INDENT}const int i_in,
                {INDENT}{DATA_TYPE} *f_out,
                {INDENT}int *i_out) {{
                    {API_NAME}(f_in,i_in,*f_out,*i_out);
                }}
                """.format(FUNC_NAME=func_name, \
                           DATA_TYPE=data_type, \
                           API_NAME=fdps_api_name, \
                           INDENT=indent_at_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_broadcast_cpp_if(self,fh):
        # Set the common indent
        indent = ""
        # Set the output texts 
        texts = []
        for key,label in self.__comm_all_types.items():
            # Get data type
            data_type = self.__ftnDT_to_cppDT[key]
            # (1) Scalar case
            func_name = "fdps_broadcast_scalar_" + label
            text = """
            void {FUNC_NAME}({DATA_TYPE} *val, int n, int src) {{
                PS::Comm::broadcast(val,n,src);
            }}
            """.format(FUNC_NAME=func_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            # (2) Array case
            func_name = "fdps_broadcast_array_" + label
            text = """
            void {FUNC_NAME}({DATA_TYPE} *val, int n, int src) {{
                PS::Comm::broadcast(val,n,src);
            }}
            """.format(FUNC_NAME=func_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_psys_fptr_meth(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the outout texts
        texts_tmp1 = []
        texts_tmp2 = []
        func_name_base = "get_psys_fptr"
        prefix_first   = "generic :: get_psys_fptr => "
        prefix_others  = " " * int(len(prefix_first))
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # (1) Private procedures
            text = "procedure, private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts_tmp1.append(text)
            # (2) Connection to generic procedure
            if (func_num == 0):
               text = prefix_first  + func_name
            else:
               text = prefix_others + func_name
            text += ", &" # continuation line mark
            text = self.__add_indent(text,indent)
            texts_tmp2.append(text) 
            func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_get_psys_fptr_decl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        func_name_base = "get_psys_fptr"
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # Set the output text
            text = "private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_psys_fptr_impl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set the output text
        texts = []
        func_name_base = "get_psys_fptr"
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name = func_name_base + func_label
            # Get the module name
            mod_name = self.__get_mod_name_by_FP(fp_t)
            # Set the output text
            text = """
            subroutine {FUNC_NAME}(this,psys_num,fptr_to_FP)
               use {MOD_NAME}
               implicit none
               class(FDPS_controller) :: this
               integer(kind=c_int), intent(IN) :: psys_num
               type({FP_T}), dimension(:), pointer, intent(INOUT) :: fptr_to_FP
               !* Local variables
               integer(kind=c_int) :: nptcl_loc
               type(c_ptr) :: cptr_to_FP

               call fdps_get_psys_cptr(psys_num,cptr_to_FP)
               nptcl_loc = fdps_get_nptcl_loc(psys_num)
               call c_f_pointer(cptr_to_FP,fptr_to_FP,[nptcl_loc])

            end subroutine {FUNC_NAME}
            """.format(FUNC_NAME=func_name,\
                       MOD_NAME=mod_name,  \
                       FP_T=fp_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_meth(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the outout texts
        texts_tmp1 = []
        texts_tmp2 = []
        func_name_base = "add_particle"
        prefix_first   = "generic :: add_particle => "
        prefix_others  = " " * int(len(prefix_first))
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # (1) Private procedures
            text = "procedure, private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts_tmp1.append(text)
            # (2) Connection to generic procedure
            if (func_num == 0):
               text = prefix_first  + func_name
            else:
               text = prefix_others + func_name
            text += ", &" # continuation line mark
            text = self.__add_indent(text,indent)
            texts_tmp2.append(text) 
            func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_add_particle_decl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        func_name_base = "add_particle"
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # Set the output text
            text = "private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_ftn_if(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        if_name_base = "fdps_add_particle"
        if_num = 0
        for fp_t in self.__FPs:
            if_name = if_name_base + "{0:05d}".format(if_num)
            mod_name = self.__get_mod_name_by_FP(fp_t)
            text = """
            subroutine {IFC_NAME}(psys_num,ptcl) bind(c)
               use, intrinsic :: iso_c_binding
               use {MOD_NAME}
               implicit none
               integer(kind=c_int), value, intent(in) :: psys_num
               type({FP_T}), intent(in) :: ptcl
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       MOD_NAME=mod_name, \
                       FP_T=fp_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            if_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_add_particle_impl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        func_name_base = "add_particle"
        func_num = 0
        for fp_t in self.__FPs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            if_name = "fdps_" + func_name
            mod_name = self.__get_mod_name_by_FP(fp_t)
            text = """
            subroutine {FUNC_NAME}(this,psys_num,ptcl)
               use, intrinsic :: iso_c_binding
               use {MOD_NAME}
               implicit none
               class(FDPS_controller) :: this
               integer(kind=c_int), intent(in) :: psys_num
               type({FP_T}), intent(in) :: ptcl

               call {IFC_NAME}(psys_num,ptcl)

            end subroutine {FUNC_NAME}
            """.format(FUNC_NAME=func_name, \
                       IFC_NAME=if_name, \
                       MOD_NAME=mod_name, \
                       FP_T=fp_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_sort_particle_ftn_if(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        if_name_base = "fdps_sort_particle"
        if_num = 0
        for fp_t in self.__FPs:
            if_name = if_name_base + "{0:05d}".format(if_num)
            mod_name = self.__get_mod_name_by_FP(fp_t)
            text = """
            subroutine {IFC_NAME}(psys_num,pfunc_comp) bind(c)
               use, intrinsic :: iso_c_binding
               use {MOD_NAME}
               implicit none
               integer(kind=c_int), value, intent(in) :: psys_num
               type(c_funptr), value, intent(in) :: pfunc_comp
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       MOD_NAME=mod_name)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            if_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_sort_particle_impl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the outout texts
        texts = []
        if_name_base = "fdps_sort_particle"
        func_num = 0
        for fp_t in self.__FPs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            if_name = if_name_base + func_label
            # Set the output text
            text = """
            case("{FP_T}")
               call {IFC_NAME}(psys_num,pfunc_comp)
            """.format(IFC_NAME=if_name, \
                       FP_T=fp_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_ftn_if(self,fh,generic_func_name):
        # Set the indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        #=================================
        # (1) Short-range interaction
        #=================================
        if_name_base = "fdps_" + generic_func_name + "_s"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}psys_num,    &
                {INDENT}dinfo_num,   &
                {INDENT}clear,       &
                {INDENT}list_mode) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep
                   logical(kind=c_bool), value, intent(in) :: clear
                   integer(kind=c_int), value, intent(in) :: list_mode
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}psys_num,    &
                {INDENT}clear) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,psys_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep
                   logical(kind=c_bool), value, intent(in) :: clear
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,search_t in self.__tree_kinds:
                if (tree_t == "Long"):
                    continue
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}dinfo_num,   &
                {INDENT}clear) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep
                   logical(kind=c_bool), value, intent(in) :: clear
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        #=================================
        # (2) Long-range interaction
        #=================================
        if_name_base = "fdps_" + generic_func_name + "_l"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}pfunc_ep_sp, &
                {INDENT}psys_num,    &
                {INDENT}dinfo_num,   &
                {INDENT}clear,       &
                {INDENT}list_mode) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,psys_num,dinfo_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
                   logical(kind=c_bool), value, intent(in) :: clear
                   integer(kind=c_int), value, intent(in) :: list_mode
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}pfunc_ep_sp, &
                {INDENT}psys_num,    &
                {INDENT}clear) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,psys_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
                   logical(kind=c_bool), value, intent(in) :: clear
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,multipole_t in self.__tree_kinds:
                if (tree_t == "Short"):
                    continue
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 12)
                # Set the output text
                text = """
                subroutine {IFC_NAME}(tree_num, &
                {INDENT}pfunc_ep_ep, &
                {INDENT}pfunc_ep_sp, &
                {INDENT}dinfo_num,   &
                {INDENT}clear) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   integer(kind=c_int), value, intent(in) :: tree_num,dinfo_num
                   type(c_funptr), value, intent(in) :: pfunc_ep_ep,pfunc_ep_sp
                   logical(kind=c_bool), value, intent(in) :: clear
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_ftn_impl_s(self,fh,generic_func_name):
        # Set the indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        if_name_base = "fdps_" + generic_func_name + "_s"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{FP_T},{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{SEARCH_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}psys_num,    &
                   {INDENT}dinfo_num,   &
                   {INDENT}clear_,      &
                   {INDENT}list_mode_)
                """.format(IFC_NAME=if_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,search_t in self.__calc_force_kinds_short:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{FP_T},{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{SEARCH_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}psys_num,    &
                   {INDENT}clear_)
                """.format(IFC_NAME=if_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_making_tree"):
            for tree_t,force_t,epi_t,epj_t,search_t in self.__tree_kinds:
                if (tree_t == "Long"):
                    continue
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{SEARCH_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}dinfo_num,   &
                   {INDENT}clear_)
                """.format(IFC_NAME=if_name, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           SEARCH_T=search_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_calc_force_ftn_impl_l(self,fh,generic_func_name):
        # Set the indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        if_name_base = "fdps_" + generic_func_name + "_l"
        func_num = 0
        if ((generic_func_name == "calc_force_all_and_write_back") or \
            (generic_func_name == "calc_force_all")):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{FP_T},{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MULTIPOLE_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}pfunc_ep_sp, &
                   {INDENT}psys_num,    &
                   {INDENT}dinfo_num,   &
                   {INDENT}clear_,      &
                   {INDENT}list_mode_)
                """.format(IFC_NAME=if_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_and_write_back"):
            for fp_t,tree_t,force_t,epi_t,epj_t,multipole_t in self.__calc_force_kinds_long:
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{FP_T},{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MULTIPOLE_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}pfunc_ep_sp, &
                   {INDENT}psys_num,    &
                   {INDENT}clear_)
                """.format(IFC_NAME=if_name, \
                           FP_T=fp_t, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        elif (generic_func_name == "calc_force_making_tree"):
            for tree_t,force_t,epi_t,epj_t,multipole_t in self.__tree_kinds:
                if (tree_t == "Short"):
                    continue
                # Set the function name
                func_label = "{0:05d}".format(func_num)
                if_name = if_name_base + func_label
                # Set the other indents
                indent_at_func_args = " " * int(len(if_name) + 6)
                # Set the output text
                text = """
                case("{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MULTIPOLE_T}")
                   call {IFC_NAME}(tree_num,    &
                   {INDENT}pfunc_ep_ep, &
                   {INDENT}pfunc_ep_sp, &
                   {INDENT}dinfo_num,   &
                   {INDENT}clear_)
                """.format(IFC_NAME=if_name, \
                           TREE_T=tree_t, \
                           FORCE_T=force_t, \
                           EPI_T=epi_t, \
                           EPJ_T=epj_t, \
                           MULTIPOLE_T=multipole_t, \
                           INDENT=indent_at_func_args)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
                func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_meth(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the output texts
        texts_tmp1 = []
        texts_tmp2 = []
        func_name_base = "get_neighbor_list"
        prefix_first   = "generic :: get_neighbor_list => "
        prefix_others  = " " * int(len(prefix_first))
        func_num = 0
        for epj_t in self.__EPJs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            # (1) Private procedures
            text = "procedure, private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts_tmp1.append(text)
            # (2) Connection to generic procedure
            if (func_num == 0):
               text = prefix_first  + func_name
            else:
               text = prefix_others + func_name
            text += ", &" # continuation line mark
            text = self.__add_indent(text,indent)
            texts_tmp2.append(text) 
            func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_get_ngb_list_decl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        func_name_base = "get_neighbor_list"
        func_num = 0
        for epj_t in self.__EPJs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # Set the output text
            text = "private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_ftn_if(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the output texts
        texts = []
        if_name_base = "fdps_get_neighbor_list"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            if (tree_t == "Long"):
                if ((mode_t != "MonopoleWithScatterSearch") and \
                    (mode_t != "QuadrupoleWithScatterSearch")):
                    continue
            # Generate
            if_name = if_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(if_name) + 12)
            text = """
            subroutine {IFC_NAME}(tree_num, &
            {INDENT}pos, &
            {INDENT}r_search, &
            {INDENT}num_epj, &
            {INDENT}cptr_to_epj) bind(c)
               use, intrinsic :: iso_c_binding
               use fdps_vector
               implicit none
               integer(kind=c_int), value, intent(in) :: tree_num
               type(fdps_f64vec), intent(in) :: pos
               real(kind=c_double), value, intent(in) :: r_search
               integer(kind=c_int), intent(inout) :: num_epj
               type(c_ptr), intent(inout) :: cptr_to_epj
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_ngb_list_impl(self,fh):
        # Set the output texts
        texts = []
        method_name_base = "get_neighbor_list"
        method_num = 0
        for target_epj_t in self.__EPJs:
            method_name = method_name_base + "{0:05d}".format(method_num)
            mod_name = self.__get_mod_name_by_EPJ(target_epj_t)
            # First-half
            indent = self.__get_ftn_indent(1)
            text = """
            subroutine {METHOD_NAME}(this,tree_num,pos,r_search,num_epj,fptr_to_EPJ)
               use {MOD_NAME}
               use fdps_vector
               implicit none
               class(FDPS_controller) :: this
               integer(kind=c_int), intent(IN) :: tree_num
               type(fdps_f64vec), intent(IN) :: pos
               real(kind=c_double), intent(IN) :: r_search
               integer(kind=c_int), intent(INOUT) :: num_epj
               type({EPJ_T}), dimension(:), pointer, intent(INOUT) :: fptr_to_EPJ
               !* Local parameters
               integer, parameter :: bufsize=256
               !* Local variables
               character(len=bufsize,kind=c_char) :: info
               type(c_ptr) :: cptr_to_EPJ
               !-(To throw errors)
               character(len=256) :: errmsg
               character(len=64) :: func_name

               call get_tree_info(this,tree_num,info)
               select case (trim(info))
            """.format(METHOD_NAME=method_name, \
                       MOD_NAME=mod_name, \
                       EPJ_T=target_epj_t)
            text = text[1:].rstrip() + "\n" 
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            # Generate case-branch
            if_name_base = "fdps_get_neighbor_list"
            if_num = 0
            for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
                # Skip if the generable condition is not satisfied.
                if (tree_t == "Long"):
                    if ((mode_t != "MonopoleWithScatterSearch") and \
                        (mode_t != "QuadrupoleWithScatterSearch")):
                        continue
                # Generate a case-branch
                if_name = if_name_base + "{0:05d}".format(if_num)
                indent_at_func_args = " " * int(len(if_name) + 6)
                if (target_epj_t == epj_t):
                    indent = self.__get_ftn_indent(2)
                    text = """
                    case("{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MODE_T}")
                       call {IFC_NAME}(tree_num, &
                       {INDENT}pos, &
                       {INDENT}r_search, &
                       {INDENT}num_epj, &
                       {INDENT}cptr_to_EPJ)
                    """.format(IFC_NAME=if_name, \
                               TREE_T=tree_t, \
                               FORCE_T=force_t, \
                               EPI_T=epi_t, \
                               EPJ_T=epj_t, \
                               MODE_T=mode_t, \
                               INDENT=indent_at_func_args)
                    text = text[1:].rstrip() + "\n" 
                    text = textwrap.dedent(text)
                    text = self.__add_indent(text,indent)
                    texts.append(text)
                if_num += 1
            # Last-half
            indent = self.__get_ftn_indent(1)
            text = """
               case default
                  errmsg =  "EssentialParticleJ specified does not have " &
                         // "a member variable representing the search "  &
                         // "radius or Tree specified does not support " &
                         // "neighbor search."
                  func_name = "{METHOD_NAME}"
                  call print_errmsg(errmsg,func_name)
                  !* Additional information
                  write(*,100)"{EPJ_T}",trim(info)
                  100 format("Please check the definitions of EssentialParticleJ"/ &
                             "and tree object:"/ &
                             " - EssentialParticleJ: ",a/ &
                             " - TreeInfo: ",a)
                  call PS_abort(this)
                  stop 1
               end select
               
               !* Convert C-pointer to Fortran-pointer
               call c_f_pointer(cptr_to_EPJ,fptr_to_EPJ,[num_epj])

            end subroutine {METHOD_NAME}
            """.format(METHOD_NAME=method_name, \
                       EPJ_T=target_epj_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            method_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_meth(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the output texts
        texts_tmp1 = []
        texts_tmp2 = []
        func_name_base = "get_epj_from_id"
        prefix_first   = "generic :: get_epj_from_id => "
        prefix_others  = " " * int(len(prefix_first))
        func_num = 0
        for epj_t in self.__EPJs:
            func_name = func_name_base + "{0:05d}".format(func_num)
            # (1) Private procedures
            text = "procedure, private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts_tmp1.append(text)
            # (2) Connection to generic procedure
            if (func_num == 0):
               text = prefix_first  + func_name
            else:
               text = prefix_others + func_name
            text += ", &" # continuation line mark
            text = self.__add_indent(text,indent)
            texts_tmp2.append(text) 
            func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_get_epj_from_id_decl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        func_name_base = "get_epj_from_id"
        func_num = 0
        for epj_t in self.__EPJs:
            # Set the function name
            func_label = "{0:05d}".format(func_num)
            func_name  = func_name_base + func_label
            # Set the output text
            text = "private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_ftn_if(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the output texts
        texts = []
        if_name_base = "fdps_get_epj_from_id"
        func_num = 0
        for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
            # Skip if the generable condition is not satisfied.
            if (not self.__getId_existence_check(epj_t)):
                continue
            # Generate
            if_name = if_name_base + "{0:05d}".format(func_num)
            indent_at_func_args = " " * int(len(if_name) + 12)
            text = """
            subroutine {IFC_NAME}(tree_num, &
            {INDENT}id, &
            {INDENT}cptr_to_epj) bind(c)
               use, intrinsic :: iso_c_binding
               use fdps_vector
               implicit none
               integer(kind=c_int), value, intent(in) :: tree_num
               integer(kind=c_long_long), value, intent(in) :: id
               type(c_ptr), intent(inout) :: cptr_to_epj
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       INDENT=indent_at_func_args)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            func_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_get_epj_from_id_impl(self,fh):
        # Set the output texts
        texts = []
        method_name_base = "get_epj_from_id"
        method_num = 0
        for target_epj_t in self.__EPJs:
            method_name = method_name_base + "{0:05d}".format(method_num)
            mod_name = self.__get_mod_name_by_EPJ(target_epj_t)
            # First-half
            indent = self.__get_ftn_indent(1)
            text = """
            subroutine {METHOD_NAME}(this,tree_num,id,fptr_to_EPJ)
               use {MOD_NAME}
               use fdps_vector
               implicit none
               class(FDPS_controller) :: this
               integer(kind=c_int), intent(IN) :: tree_num
               integer(kind=c_long_long), intent(IN) :: id
               type({EPJ_T}), pointer, intent(INOUT) :: fptr_to_EPJ
               !* Local parameters
               integer, parameter :: bufsize=256
               !* Local variables
               character(len=bufsize,kind=c_char) :: info
               type(c_ptr) :: cptr_to_EPJ
               !-(To throw errors)
               character(len=256) :: errmsg
               character(len=64) :: func_name

               call get_tree_info(this,tree_num,info)
               select case (trim(info))
            """.format(METHOD_NAME=method_name, \
                       MOD_NAME=mod_name, \
                       EPJ_T=target_epj_t)
            text = text[1:].rstrip() + "\n" 
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            # Generate case-branch
            if_name_base = "fdps_get_epj_from_id"
            if_num = 0
            for tree_t,force_t,epi_t,epj_t,mode_t in self.__tree_kinds:
                # Skip if the generable condition is not satisfied.
                if (not self.__getId_existence_check(epj_t)):
                    continue
                # Generate a case-branch
                if_name = if_name_base + "{0:05d}".format(if_num)
                indent_at_func_args = " " * int(len(if_name) + 6)
                if (target_epj_t == epj_t):
                    indent = self.__get_ftn_indent(2)
                    text = """
                    case("{TREE_T},{FORCE_T},{EPI_T},{EPJ_T},{MODE_T}")
                       call {IFC_NAME}(tree_num, &
                       {INDENT}id, &
                       {INDENT}cptr_to_EPJ)
                    """.format(IFC_NAME=if_name, \
                               TREE_T=tree_t, \
                               FORCE_T=force_t, \
                               EPI_T=epi_t, \
                               EPJ_T=epj_t, \
                               MODE_T=mode_t, \
                               INDENT=indent_at_func_args)
                    text = text[1:].rstrip() + "\n" 
                    text = textwrap.dedent(text)
                    text = self.__add_indent(text,indent)
                    texts.append(text)
                if_num += 1
            # Last-half
            indent = self.__get_ftn_indent(1)
            text = """
               case default
                  errmsg =  "EssentialParticleJ specified does not have " &
                         // "a method returning the particle id, or, "  &
                         // "unknown tree object is specified."
                  func_name = "{METHOD_NAME}"
                  call print_errmsg(errmsg,func_name)
                  !* Additional information
                  write(*,100)"{EPJ_T}",trim(info)
                  100 format("Please check the definitions of EssentialParticleJ"/ &
                             "and tree object:"/ &
                             " - EssentialParticleJ: ",a/ &
                             " - TreeInfo: ",a)
                  call PS_abort(this)
                  stop 1
               end select
               
               !* Convert C-pointer to Fortran-pointer
               call c_f_pointer(cptr_to_EPJ,fptr_to_EPJ)

            end subroutine {METHOD_NAME}
            """.format(METHOD_NAME=method_name, \
                       EPJ_T=target_epj_t)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            method_num += 1
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_reduction_meth(self,fh,generic_func_name,gen_func_w_id):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Make the list of function names
        func_names = []
        for label in self.__comm_all_types.values():
            func_name  = generic_func_name + "_" + label
            func_names.append(func_name)
        if (gen_func_w_id):
            for label in self.__comm_real_types.values():
                func_name  = generic_func_name + "_w_id_" + label
                func_names.append(func_name)
        # Set the output text
        texts_tmp1 = []
        texts_tmp2 = []
        prefix_first  = "generic :: " + generic_func_name + " => "
        prefix_others = " " * int(len(prefix_first))
        func_num = 0
        for func_name in func_names:
            # Private procedures 
            text = "procedure, private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts_tmp1.append(text)
            # Connection to generic procedure
            if (func_num == 0):
               text = prefix_first  + func_name
            else:
               text = prefix_others + func_name
            text += ", &" # continuation line mark
            text = self.__add_indent(text,indent)
            texts_tmp2.append(text) 
            func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_reduction_decl(self,fh,generic_func_name,gen_func_w_id):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Make the list of function names
        func_names = []
        for label in self.__comm_all_types.values():
            func_name = generic_func_name + "_" + label
            func_names.append(func_name)
        if (gen_func_w_id):
            for label in self.__comm_real_types.values():
                func_name = generic_func_name + "_w_id_" + label
                func_names.append(func_name)
        # Set the outout texts
        texts = []
        for func_name in func_names:
            text = "private :: {0}".format(func_name)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_reduction_ftn_if(self,fh,generic_func_name,gen_func_w_id):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        # (1) Case: w/o ID
        for data_type,label in self.__comm_all_types.items():
            if_name = "fdps_" + generic_func_name + "_" + label
            text = """
            subroutine {IFC_NAME}(f_in,f_out) bind(c)
               use, intrinsic :: iso_c_binding
               implicit none
               {DATA_TYPE}, value, intent(in) :: f_in
               {DATA_TYPE}, intent(inout) :: f_out
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) Case: w/ ID
        if (gen_func_w_id):
            for data_type,label in self.__comm_real_types.items():
                if_name = "fdps_" + generic_func_name + "_w_id_" + label
                text = """
                subroutine {IFC_NAME}(f_in,i_in,f_out,i_out) bind(c)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   {DATA_TYPE}, value, intent(in) :: f_in
                   {DATA_TYPE}, intent(inout) :: f_out
                   integer(kind=c_int), value, intent(in) :: i_in
                   integer(kind=c_int), intent(inout) :: i_out
                end subroutine {IFC_NAME}
                """.format(IFC_NAME=if_name, \
                           DATA_TYPE=data_type)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_reduction_ftn_impl(self,fh,generic_func_name,gen_func_w_id):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Case: w/o ID
        for data_type,label in self.__comm_all_types.items():
            func_name = generic_func_name + "_" + label
            if_name = "fdps_" + generic_func_name + "_" + label
            text = """
            subroutine {FUNC_NAME}(this,f_in,f_out)
               use, intrinsic :: iso_c_binding
               implicit none
               class(FDPS_controller) :: this
               {DATA_TYPE}, intent(IN) :: f_in
               {DATA_TYPE}, intent(INOUT) :: f_out

                call {IFC_NAME}(f_in,f_out)

            end subroutine {FUNC_NAME}
            """.format(FUNC_NAME=func_name, \
                       IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # (2) Case: w/ ID
        if (gen_func_w_id):
            for data_type,label in self.__comm_real_types.items():
                func_name = generic_func_name + "_w_id_" + label
                if_name = "fdps_" + generic_func_name + "_w_id_" + label
                text = """
                subroutine {FUNC_NAME}(this,f_in,i_in,f_out,i_out)
                   use, intrinsic :: iso_c_binding
                   implicit none
                   class(FDPS_controller) :: this
                   {DATA_TYPE}, intent(IN) :: f_in
                   {DATA_TYPE}, intent(INOUT) :: f_out
                   integer(kind=c_int), intent(IN) :: i_in
                   integer(kind=c_int), intent(INOUT) :: i_out

                   call {IFC_NAME}(f_in,i_in,f_out,i_out)

                end subroutine {FUNC_NAME}
                """.format(FUNC_NAME=func_name, \
                           IFC_NAME=if_name, \
                           DATA_TYPE=data_type)
                text = text[1:].rstrip() + "\n"
                text = textwrap.dedent(text)
                text = self.__add_indent(text,indent)
                texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_broadcast_meth(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set the outout texts
        texts_tmp1 = []
        texts_tmp2 = []
        prefix_first   = "generic :: broadcast => "
        prefix_others  = " " * int(len(prefix_first))
        func_num = 0
        for label in self.__comm_all_types.values():
            # Make the list of function names
            func_names = []
            func_names.append("broadcast_scalar_" + label)
            func_names.append("broadcast_array_"  + label)
            for func_name in func_names:
                # Private procedures 
                text = "procedure, private :: {0}".format(func_name)
                text = self.__add_indent(text,indent)
                texts_tmp1.append(text)
                # Connection to generic procedure
                if (func_num == 0):
                   text = prefix_first  + func_name
                else:
                   text = prefix_others + func_name
                text += ", &" # continuation line mark
                text = self.__add_indent(text,indent)
                texts_tmp2.append(text) 
                func_num += 1
        # Output the text
        text = "".join(texts_tmp1 + texts_tmp2)
        end = text.rindex(", &")
        text = text[:end] + "\n"
        fh.write(text)

    def __write_broadcast_decl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        func_num = 0
        for label in self.__comm_all_types.values():
            # Make the list of function names
            func_names = []
            func_names.append("broadcast_scalar_" + label)
            func_names.append("broadcast_array_" + label)
            for func_name in func_names:
                # Set the output text
                text = "private :: {0}".format(func_name)
                text = self.__add_indent(text,indent)
                texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_broadcast_ftn_if(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(2)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        for data_type,label in self.__comm_all_types.items():
            # (1) Case: scalar argument
            if_name = "fdps_broadcast_scalar_" + label
            text = """
            subroutine {IFC_NAME}(val,n,src) bind(c)
               use, intrinsic :: iso_c_binding
               implicit none
               {DATA_TYPE}, intent(inout) :: val
               integer(kind=c_int), value, intent(in) :: n,src
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            # (2) Case: array argument
            if_name = "fdps_broadcast_array_" + label
            text = """
            subroutine {IFC_NAME}(vals,n,src) bind(c)
               use, intrinsic :: iso_c_binding
               implicit none
               integer(kind=c_int), value, intent(in) :: n
               {DATA_TYPE}, dimension(n), intent(inout) :: vals
               integer(kind=c_int), value, intent(in) :: src
            end subroutine {IFC_NAME}
            """.format(IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def __write_broadcast_ftn_impl(self,fh):
        # Set the common indent
        indent = self.__get_ftn_indent(1)
        # Set th outout texts
        texts = []
        # (1) Private procedures
        for data_type,label in self.__comm_all_types.items():
            # (1) Scalar case
            func_name = "broadcast_scalar_" + label
            if_name = "fdps_broadcast_scalar_" + label
            text = """
            subroutine {FUNC_NAME}(this,val,n,src)
               use, intrinsic :: iso_c_binding
               implicit none
               class(FDPS_controller) :: this
               {DATA_TYPE}, intent(INOUT) :: val
               integer(kind=c_int), intent(in) :: n
               integer(kind=c_int), optional, intent(IN) :: src

                if (present(src)) then
                   call {IFC_NAME}(val,1,src)
                else
                   call {IFC_NAME}(val,1,0)
                end if

            end subroutine {FUNC_NAME}
            """.format(FUNC_NAME=func_name, \
                       IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
            # (2) Array case
            func_name = "broadcast_array_" + label
            if_name = "fdps_broadcast_array_" + label
            text = """
            subroutine {FUNC_NAME}(this,vals,n,src)
               use, intrinsic :: iso_c_binding
               implicit none
               class(FDPS_controller) :: this
               integer(kind=c_int), intent(in) :: n
               {DATA_TYPE}, dimension(n), intent(INOUT) :: vals
               integer(kind=c_int), optional, intent(IN) :: src
            
               if (present(src)) then
                  call {IFC_NAME}(vals,n,src)
               else
                  call {IFC_NAME}(vals,n,0)
               end if

            end subroutine {FUNC_NAME}
            """.format(FUNC_NAME=func_name, \
                       IFC_NAME=if_name, \
                       DATA_TYPE=data_type)
            text = text[1:].rstrip() + "\n"
            text = textwrap.dedent(text)
            text = self.__add_indent(text,indent)
            texts.append(text)
        # Output the text
        text = "".join(texts)
        fh.write(text)

    def read_files(self,file_names):
        # Output to stdout
        self.__print_checkpoint("reading user's Fortran files")
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
            # Compile regexp pattern
            pat_ftn_comment = re.compile(r"^!.*")
            pat_fdps_dir    = re.compile(r"^!\$fdps")
            # Read the file
            buf_addr = 0
            line_prev = ""
            for line in fp:
                # Delete LF
                line = line.rstrip()
                # Delete the leading and trailing blank characters in the line
                line = line.strip() 
                # Delete tab characetrs
                line = line.strip("\t")
                # Skip if the line contains nothing
                if (line == ""): continue
                # Skip if the line contains Fortran comment only
                if (pat_ftn_comment.match(line)):
                    if (pat_fdps_dir.match(line) is None):
                        continue
                # Delete the continuation line mark (&) at the head of the line
                # [*** Note ***]
                #    Because we already delete the leading blank characters,
                #    `&` mark must exist at head of the line if exists.
                Mobj = re.match(r"^&(?P<line_rmng>.*)",line)
                if (Mobj):
                    line = Mobj.group('line_rmng')
                    line = line.strip()
                # Connect continuation line to the previous line
                if (line_prev is not ""):
                    if (pat_ftn_comment.match(line)):
                        # This is invalid beacause we cannot connect a line
                        # that consists of Fortran comment or FDPS directive
                        # to the previous one.
                        self.__print_error("An invalid & mark is detected at the line:")
                        print("{0}".format(line))
                        sys.exit()
                    else:
                        # In this case, the current line is connected to
                        # the previous line.
                        line = line_prev + line
                # Detect the continuation line mark (&) at the end of the line
                Mobj = re.match(r"^(?P<line_body>[^!]*)&",line)
                if (Mobj):
                    # In this case, this line has a continuation line
                    line = Mobj.group('line_body')
                    line = line.strip()
                    line_prev = line
                    continue
                else:
                    # In this case, the line does not have a continuation line.
                    # We reset line_prev to empty.
                    line_prev = ""
                # Save the line 
                f.content.append([buf_addr,line])
                buf_addr += 1
            # Close the file
            fp.close()

    def analyze(self):
        # Output to stdout
        self.__print_checkpoint("analyzing Fortran files ...")
        # Compile regexp patterns
        pat_module_start  = re.compile(r"(?:^module.*)\Z")
        pat_module_end    = re.compile(r"(?:^end module.*)\Z")
        pat_struct_start  = re.compile(r"(?:^type.*)\Z")
        pat_struct_end    = re.compile(r"(?:^end type.*)\Z")
        pat_struct_target = re.compile(r".*!\$fdps.*(FP|EPI|EPJ|Force)")
        # Loop w.r.t. files
        print("identifying modules and structures ...")
        for f in self.files:
            # Check the file name
            print("processing ... the file {0}".format(f.name))
            # Initialize flags and stacks
            module_detect_mode = False
            struct_detect_mode = False
            module_stack = []
            struct_stack = []
            #=============================================================================
            # Sweep the file content to detect modules and structures
            line_prev = ""
            for item in f.content:
                # Extract informaiton
                buf_addr = item[0]
                line     = item[1]
                # Seek keywords `module` or `type`
                if (module_detect_mode == False):
                    if (pat_module_start.match(line)):
                        # Get the module name
                        module_name = re.match(r"^module (?P<module_name>.*)",line).group('module_name')
                        # Create a module object
                        f.modules.append(Module_Data())
                        # Register to the stack
                        module_stack.append([buf_addr,module_name])
                        # Rise module_detect_mode
                        module_detect_mode = True
                if (module_detect_mode == True):
                    if (pat_module_end.match(line)):
                        # Get the reference for the current module
                        mod = f.modules[-1]
                        # Setup the module object
                        mod.buf_addr_sta, mod.name = module_stack.pop()
                        mod.buf_addr_end = buf_addr
                        # Lower the flag
                        module_detect_mode = False
                if (struct_detect_mode == False):
                    if (pat_struct_start.match(line)):
                        if (pat_struct_target.match(line) or \
                            pat_struct_target.match(line_prev)):
                            # Get the structure name
                            struct_name = re.match(r"type.*::\s*(?P<struct_name>\w+)",line).group('struct_name')
                            struct_name = struct_name.lower() # Enforce to lower-case
                            # Get the reference for the current module
                            mod = f.modules[-1]
                            # Create a struct object
                            mod.structures.append(Structure_Data())
                            s = mod.structures[-1]
                            # Setup attributes
                            for key in self.__usable_fdps_str_dirs:
                                if (re.search(key,line) or \
                                    re.search(key,line_prev)):
                                    s.attrib[key] = True
                            # Register to the stack
                            struct_stack.append([buf_addr,struct_name])
                            # Raise the flag
                            struct_detect_mode = True
                if (struct_detect_mode == True):
                    if (pat_struct_end.match(line)):
                        # Get the reference for the current module
                        mod = f.modules[-1]
                        # Get the reference for the current structure
                        s = mod.structures[-1]
                        # Setup struct
                        s.buf_addr_sta, s.name = struct_stack.pop()
                        s.buf_addr_end = buf_addr
                        # Lower the flag
                        struct_detect_mode = False
                # Set the previous line
                line_prev = line

            #=============================================================================
            # Identify member variables
            print("identifying member variables in each structure ...")
            pat_fdps_dir_only = re.compile(r"(?:^!\$fdps.*)\Z")
            for mod in f.modules:
                for s in mod.structures:
                    buf_addr_sta = s.buf_addr_sta + 1
                    buf_addr_end = s.buf_addr_end - 1
                    for buf_addr in range(buf_addr_sta,buf_addr_end+1):
                        line      = (f.content[buf_addr])[1]
                        line_prev = (f.content[buf_addr-1])[1]
                        # Skip if the line consists of FDPS directives only
                        if (pat_fdps_dir_only.match(line)):
                            continue
                        # Analyze the member variable
                        #-----------------------------------------------------------------
                        # (1) Examine if the line has the `::` symbol or not.
                        #     If exits, we divide the line by the symbol.
                        #     Otherwise, we divide the line according to classical
                        #     Fortran rule.
                        if (re.match(r"[^:]*::[^:]*",line)):
                            # In this case, the `::` symbol exists.
                            Mobj = re.match(r"^(?P<lhs>.*)::(?P<rhs>.*)",line)
                        else:
                            # In this case, the `::` symbol does not exist.
                            Mobj = re.match(r"^(?P<lhs>[a-zA-Z0-9_()=,]*\))\s+(?P<rhs>[a-zA-Z].*)",line)
                            # [*** Note ***]
                            #    When using iso_c_binding module, all the data types in a 
                            #    structure must be interoperable with C.
                            #    Therefore, a contiguious space that separates lhs from rhs
                            #    is bounded by the symbol `)` at its left side. 
                            #    (note that `dimension(*)` also ends with `)`.)
                            #    The right side of the space(s) must be bounded by 
                            #    a character in [a-zA-Z] because rhs lists Fortran variables.
                        lhs = Mobj.group('lhs')
                        lhs = lhs.strip()
                        rhs = Mobj.group('rhs')
                        rhs = rhs.strip()
                        if ((lhs == "") or (rhs == "")):
                            self.__print_error("Syntax error detected at the line:")
                            print("{0}".format(line))
                            sys.exit()
                        # [Option] Check
                        #print("lhs --- {0}".format(lhs))
                        #print("rhs --- {0}".format(rhs))
                        #-----------------------------------------------------------------
                        # (2) Analyze rhs
                        #    (i)   Examine # of member variables.
                        #    (ii)  Check if each member is a valid array or not.
                        #    (iii) Check if the description of FDPS directive is correct
                        #          if it exists.
                        Mobj = re.search(r"!\$fdps",rhs)
                        if (Mobj is None):
                            # In this case, there is no FDPS directive in rhs.
                            mbr_defs = rhs.split(',')
                            fdps_dir = ""
                            # [*** Note ***]
                            #    We assume here that rhs does not contain multi-dimensional
                            #    Fortran arrays because it is not self-evident how to convert
                            #    them into C++ data. In this case, spliting by comma will work.
                        else:
                            # In this case, there is FDPS directive in rhs
                            Mobj = re.match(r"^(?P<mbr_defs>.*)!\$fdps(?P<fdps_dir>.*)",rhs)
                            mbr_defs = Mobj.group('mbr_defs')
                            mbr_defs = mbr_defs.strip()
                            mbr_defs = mbr_defs.split(',')
                            fdps_dir = Mobj.group('fdps_dir')
                            fdps_dir = fdps_dir.strip().lower()
                        # Examine # of member variables and check the syntax
                        members_cand = []
                        for item in mbr_defs:
                            # Delete extra spaces
                            item = item.strip()
                            # Create a tentative member object
                            members_cand.append(Member_Data())
                            mbr = members_cand[-1]
                            # Check if it is an array 
                            Mobj = re.search(r"\((?P<array_dim>.+)\)",item)
                            if (Mobj is not None):
                                is_array  = True
                                array_dim = Mobj.group('array_dim')
                                if (array_dim.isdigit() == False):
                                    self.__print_error("an invalid array description: {0}".format(item))
                                    sys.exit()
                            else:
                                is_array  = False
                                array_dim = -1
                            # Check if it has correct name
                            if (is_array):
                                Mobj = re.match(r"(?P<name>\w+)",item)
                                name = Mobj.group('name')
                            else:
                                name = item
                            name = name.lower() # Enforce to lower-case
                            if (re.match(r"^[a-zA-Z0-9_]+$",name) is None):
                                self.__print_error("incorrect name : {0}".format(item))
                                sys.exit()
                            # Save the data
                            mbr.name = name
                            mbr.is_array = is_array
                            mbr.array_dim = array_dim
                        # Check the syntax of FDPS directive if it exists 
                        if (fdps_dir != ""):
                            # Cannot associate FDPS directive with members more than one
                            if (len(members_cand) > 1):
                                self.__print_error("FDPS directive can modify a single member!")
                                sys.exit()
                            if (fdps_dir in self.__usable_fdps_mbr_dirs):
                                attrib = fdps_dir
                            else:
                                self.__print_error("unknown directive: {0}".format(fdps_dir))
                                sys.exit()
                        else:
                            attrib = ""
                        # Save the attribute
                        for mbr in members_cand:
                            mbr.attrib = attrib
                        #-----------------------------------------------------------------
                        # (3) Analyze lhs
                        #    (i)  Extract data type and check it
                        #    (ii) Check if there is dimension(*) statement and check if
                        #         it is described in the correct syntax
                        # Examine data type
                        items = lhs.split(",")
                        if (items[0].lower() not in self.__usable_data_types):
                            msg = "An invalid data type `{0}` is used in your Fortran code!".format(items[0])
                            self.__print_error(msg)
                            sys.exit()
                        else:
                            data_type = items[0].lower()
                        # Check the description of `dimension(*)` 
                        Mobj = re.search(r"dimension\((?P<array_dim>.+)\)",lhs,re.IGNORECASE)
                        if (Mobj is not None):
                            # In this case, this member is an array
                            is_array = True
                            array_dim = Mobj.group('array_dim')
                            if (array_dim.isdigit()):
                                # In this case, the dimension of the array is given by a number
                                mbr.array_dim = int(array_dim)
                            else:
                                self.__print_error("The dimension of an array is not a number!")
                                sys.exit()
                            # Check if rhs contains an array (this is a error case)
                            for mbr in members_cand:
                                if (mbr.is_array):
                                    msg = "cannot use `dimension(*)` and `var_name(*)` " \
                                        + "simultaneously to declare an array."
                                    self.__print_error(msg)
                                    sys.exit()
                        else:
                            is_array  = False
                            array_dim = -1
                        # Save the result
                        for mbr in members_cand:
                            mbr.buf_addr  = buf_addr
                            mbr.data_type = data_type
                            if (is_array):
                                mbr.is_array  = True
                                mbr.array_dim = int(array_dim)
                        #-----------------------------------------------------------------
                        # (4) Check if the previous line contains FDPS directive
                        if (pat_fdps_dir_only.match(line_prev)):
                            Mobj = re.match(r"^!\$fdps\s+(?P<key>\w+)",line_prev)
                            key  = Mobj.group('key')
                            if (key in self.__usable_fdps_mbr_dirs):
                                # Produce an error when 
                                if (len(members_cand) > 1):
                                    self.__print_error("FDPS directive can modify a single member!")
                                # Check if mbr.attrib is already defined 
                                mbr = members_cand[-1] # because of sigle member
                                if (mbr.attrib == ""):
                                    mbr.attrib = key
                                else:
                                    if (mbr.attrib is not key):
                                        self.__print_error("find inconsistent directives")
                                        sys.exit()
                            else:
                                if (key not in self.__usable_fdps_meth_dirs):
                                    self.__print_warning("ignore unknown directive:{0}".format(key))
                        #-----------------------------------------------------------------
                        # (5) Register to the current structure
                        for mbr in members_cand:
                            s.members.append(mbr)
                            #print("name      : " + mbr.name)
                            #print("attrib    : " + mbr.attrib)
                            #print("buf_addr  : " + str(mbr.buf_addr))
                            #print("data_type : " + mbr.data_type)
                            #print("is_array  : " + mbr.is_array)
                            #print("array_dim : " + mbr.array_dim)
                            #print("^^^")
                        #-----------------------------------------------------------------

            #=============================================================================
            # Check the consistency of member variables
            #   (i)  Check the data types of position,velocity,charge,rsearch
            #   (ii) Check if position member is unique
            print("checking the consistency of members variables")
            for mod in f.modules:
                for s in mod.structures:
                    fdps_members = {"id":False, \
                                    "position":False, \
                                    "velocity":False, \
                                    "charge":False, \
                                    "rsearch":False}
                    for mbr in s.members:
                        if (mbr.attrib == "id"):
                            # Check if the FDPS member directive is multiply defined
                            if (fdps_members[mbr.attrib] == False):
                                fdps_members[mbr.attrib] == True
                            else:
                                msg = "The FDPS directive `{0}` multiply defined.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()
                            # Check if data type is correct
                            if (not self.__is_s64(mbr)):
                                msg = "The data type of the member `{0}` ".format(mbr.name) \
                                    + "of the structure `{0}` ".format(s.name) \
                                    + "is not correct for `{0}`.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()
                        elif ((mbr.attrib == "position") or \
                              (mbr.attrib == "velocity")):
                            # Check if the FDPS member directive is multiply defined
                            if (fdps_members[mbr.attrib] == False):
                                fdps_members[mbr.attrib] == True
                            else:
                                msg = "The FDPS directive `{0}` multiply defined.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()
                            # Check if data type is correct
                            if (not (self.__is_f32vec(mbr) or self.__is_f64vec(mbr))):
                                msg = "The data type of the member `{0}` ".format(mbr.name) \
                                    + "of the structure `{0}` ".format(s.name) \
                                    + "is not correct vector for `{0}`.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()
                        elif ((mbr.attrib == "charge") or \
                              (mbr.attrib == "rsearch")):
                            # Check if the FDPS member directive is multiply defined
                            if (fdps_members[mbr.attrib] == False):
                                fdps_members[mbr.attrib] == True
                            else:
                                msg = "The FDPS directive `{0}` multiply defined.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()
                            # Check if data type is correct
                            if (not (self.__is_f32(mbr) or self.__is_f64(mbr))):
                                msg = "The data type of the member `{0}` ".format(mbr.name) \
                                    + "of the structure `{0}` ".format(s.name) \
                                    + "is incorrect for `{0}`.".format(mbr.attrib)
                                self.__print_error(msg)
                                sys.exit()

            #=============================================================================
            # Make the list of generable methods (s.method_DB)
            print("making method_DB")
            pat_fdps_meth_dir = re.compile(r"(?:^!\$fdps\s+(copyFromForce|copyFromForcePM|copyFromFP|clear).*)\Z")
            for mod in f.modules:
                for s in mod.structures:
                    # (1) Check FDPS member directives
                    for mbr in s.members:
                        data_type, array_dim = self.__get_fdpsDT(mbr)
                        if (mbr.attrib == "id"):
                            (s.method_DB["getId"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["getId"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                        elif (mbr.attrib == "position"):
                            (s.method_DB["getPos"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["getPos"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                            (s.method_DB["setPos"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["setPos"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                        elif (mbr.attrib == "charge"):
                            (s.method_DB["getCharge"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["getCharge"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                            (s.method_DB["getChargePM"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["getChargePM"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                        elif (mbr.attrib == "rsearch"):
                            (s.method_DB["getRSearch"])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB["getRSearch"])[s._INDX_DATA_TO_GEN]    = [data_type,mbr.name]
                    # (2) Check FDPS method directives
                    buf_addr_sta = s.buf_addr_sta + 1
                    buf_addr_end = s.buf_addr_end - 1
                    for buf_addr in range(buf_addr_sta,buf_addr_end+1):
                        line = (f.content[buf_addr])[1]
                        # Check if the line has FDPS method directives
                        if (pat_fdps_meth_dir.match(line)):
                            # Extract the class name and member pairs
                            items = line.split()
                            # Check the directive keywords
                            keyword = items[1]
                            if (keyword in self.__usable_fdps_meth_dirs):
                                if ((keyword == "copyFromForce") or \
                                    (keyword == "copyFromFP")):
                                    # Check the # of items (at least 4 items must exist)
                                    num_items = len(items)
                                    if (num_items < 4):
                                        self.__print_error("lack of information: {0}".format(keyword))
                                        sys.exit()
                                    # Get the class name
                                    class_name = items[2]
                                    class_name = class_name.lower() # Enforce to lower-case
                                    # Get the member pairs
                                    listtmp = []
                                    for Mobj in re.finditer(r"\((?P<dst_name>\w+)\s*,\s*(?P<src_name>\w+)\)",line):
                                        dst_name = Mobj.group('dst_name')
                                        src_name = Mobj.group('src_name')
                                        if ((dst_name == "") or (src_name == "")):
                                            msg = "syntax error in {0}".format(keyword)
                                            self.__print_error(msg)
                                            sys.exit()
                                        # Enforce to lower-case
                                        dst_name = dst_name.lower()
                                        src_name = src_name.lower()
                                        # Check if the member `dst_name` really exist?
                                        is_exist = False
                                        for mbr in s.members:
                                            if (mbr.name == dst_name):
                                                is_exist = True
                                        if (is_exist == False):
                                            self.__print_error("the member {0} does not exist.".format(dst_name))
                                            sys.exit()
                                        # Save the result 
                                        listtmp.append([dst_name,src_name])
                                    # Set the result
                                    obj = Method_Generation_Data()
                                    obj.class_name   = class_name
                                    obj.member_pairs = listtmp
                                    (s.method_DB[keyword])[s._INDX_GENERABLE_FLAG] = True
                                    (s.method_DB[keyword])[s._INDX_DATA_TO_GEN].append(obj)
                                elif (keyword == "copyFromForcePM"):
                                    # Get the member name
                                    mbr_name = items[2]
                                    mbr_name = mbr_name.lower()
                                    # Check if it really exists or not
                                    is_exist = False
                                    for mbr in s.members:
                                        if (mbr.name == mbr_name):
                                            is_exist = True
                                    if (is_exist == False):
                                        self.__print_error("the member {0} does not exist".format(mbr_name))
                                        sys.exit()
                                    # Save the result
                                    (s.method_DB["copyFromForcePM"])[s._INDX_GENERABLE_FLAG] = True
                                    (s.method_DB["copyFromForcePM"])[s._INDX_DATA_TO_GEN]    = mbr_name
                                elif (keyword == "clear"):
                                    # Examine which format is used
                                    num_items = len(items)
                                    if ((num_items >= 3) and (items[2] == "subroutine")):
                                        # In this case, the clear method calls the Fortran subroutine
                                        # specified by an user.
                                        gen_mode = 1
                                        # (1) Get the name of Fortran subroutine
                                        if (num_items > 3):
                                            subroutine_name = items[3]
                                            if (re.match(r"[0-9a-zA-Z_]+",subroutine_name) is None):
                                                msg = "The name of Fortran subroutine violates Fortran syntax!"
                                                self.__print_error(msg)
                                                sys.exit()
                                            # Enforce to lower-case
                                            subroutine_name = subroutine_name.lower()
                                        else:
                                            msg = "The name of Fortran subroutine is not given!"
                                            self.__print_error(msg)
                                            sys.exit()
                                        # (2) Set method_DB
                                        (s.method_DB[keyword])[s._INDX_GENERABLE_FLAG] = True
                                        (s.method_DB[keyword])[s._INDX_DATA_TO_GEN] = [gen_mode,subroutine_name]
                                    else:
                                        # In this case, the default initial values are used in the clear method
                                        # except for the member variables described after the `clear` keyword.
                                        gen_mode = 0
                                        # (1) First, we apply the default initialization to all the member
                                        #     variables of this derived data type.
                                        members = []
                                        data_types = []
                                        array_dims = []
                                        vals = []
                                        for mbr in s.members:
                                            fdpsDT, array_dim = self.__get_fdpsDT(mbr)
                                            val = self.__default_IV_of_cppDT[fdpsDT]
                                            members.append(mbr.name)
                                            data_types.append(fdpsDT)
                                            array_dims.append(array_dim)
                                            vals.append(val)
                                        # (2) Next, we take into account the exceptions.
                                        for Mobj in re.finditer(r"""(?P<mbr_name>\w+)\s*=\s*(?P<val>([^[\](){}!"'#$%&=~^|\`@*?/<>,:;]+|"([^"]|(""))+"|'([^']|(''))+'))\s*(,|$)""",line):
                                            # (2-1) Extract `mbr` and `val`
                                            mbr_name = Mobj.group('mbr_name')
                                            val_ftn  = Mobj.group('val')
                                            #print("mbr_name = {0}".format(mbr_name))
                                            #print("val_ftn  = {0}".format(val_ftn))
                                            # (2-2) Check the basic syntax
                                            if ((mbr_name == "") or (val_ftn == "")):
                                                msg = "syntax error in FDPS directive for clear!"
                                                self.__print_error(msg)
                                                sys.exit()
                                            mbr_name = mbr_name.lower() # Enforce to lower-case
                                            # (2-3) Check if mbr_name exists or not.
                                            is_exist = False
                                            for mbr in s.members:
                                                if (mbr_name == mbr.name):
                                                   is_exist = True
                                            if (is_exist == False):
                                                   msg = "{0} is not a member variable!".format(mbr_name)
                                                   self.__print_error(msg)
                                                   sys.exit()
                                            # (2-4) Check if val is keep or not
                                            if (val_ftn == "keep"):
                                                # In this case, we must skip the initialization of
                                                # this member variables. Hence, we delete it from
                                                # members, array_dims, and vals.
                                                for i in range(0,len(members)):
                                                    if (mbr_name == members[i]):
                                                        del members[i]
                                                        del data_types[i]
                                                        del array_dims[i]
                                                        del vals[i]
                                                        break
                                            else:
                                                # In this case, we update the initial value of this
                                                # member variables. In order to do it, we convert
                                                # Fortran values to C++ values.
                                                if (re.match(r"^([+]|[-])?[0-9]+$",val_ftn)):
                                                    # Integer value
                                                    val = val_ftn
                                                elif (re.match(r"^([+]|[-])?[0-9]+[.]?[0-9]*(e([+]|[-])?[0-9]+)?$",val_ftn,flags=re.IGNORECASE)):
                                                    # Floating point value (32bit)
                                                    # [ex.]
                                                    #    0.
                                                    #    1.0
                                                    #    1.0e0
                                                    val = val_ftn
                                                elif (re.match(r"^([+]|[-])?[0-9]+[.]?[0-9]+d([+]|[-])?[0-9]+$",val_ftn,flags=re.IGNORECASE)):
                                                    # Floating point value (64bit)
                                                    # [ex.]
                                                    #    1.0d0
                                                    val = re.sub(r"(d|D)","e",val_ftn)
                                                elif (re.match(r"^[.]true[.]$",val_ftn,flags=re.IGNORECASE)):
                                                    # Logical value (.true.)
                                                    val = re.sub(r"[.]true[.]","true",val_ftn,flags=re.IGNORECASE)
                                                elif (re.match(r"^[.]false[.]$",val_ftn,re.IGNORECASE)):
                                                    # Logical value (.false.)
                                                    val = re.sub(r"[.]false[.]","false",val_ftn,flags=re.IGNORECASE)
                                                elif (re.match(r"^'.*'$",val_ftn)):
                                                    # Fortran character constant quoted by single quotes
                                                    Mobj = re.match(r"^'(?P<val>.*)'$",val_ftn)
                                                    val = Mobj.group('val')
                                                    val = re.sub(r"\\",r"\\\\",val)
                                                    val = re.sub(r"\?",r"\?",val)
                                                    val = re.sub(r'"',r'\"',val)
                                                    val = re.sub(r"''",r"\'",val)
                                                    val = '"' + val + '"'
                                                elif (re.match(r'^".*"$',val_ftn)):
                                                    # Fortran character constant quoted by double quotes
                                                    Mobj = re.match(r'^"(?P<val>.*)"$',val_ftn)
                                                    val = Mobj.group('val')
                                                    val = re.sub(r"\\",r"\\\\",val)
                                                    val = re.sub(r"\?",r"\?",val)
                                                    val = re.sub(r"'",r"\'",val)
                                                    val = re.sub(r'""',r'\"',val)
                                                    val = '"' + val + '"'
                                                else:
                                                    # Unsupported cases
                                                    msg = "Unknown value is specified in FDPS directive clear!"
                                                    self.__print_error(msg)
                                                    sys.exit()
                                                # Update
                                                for i in range(0,len(members)):
                                                    if (mbr_name == members[i]):
                                                        vals[i] = val
                                        # (3) Set method_DB
                                        (s.method_DB[keyword])[s._INDX_GENERABLE_FLAG] = True
                                        (s.method_DB[keyword])[s._INDX_DATA_TO_GEN] = [gen_mode,members,data_types,array_dims,vals]
                            else:
                                self.__print_warning("ignore an unknown FDPS directive")
                    # (3) Detect the case that the description !$fdps clear does not exist.
                    if (s.attrib["Force"] == True):
                        keyword = "clear"
                        data_to_gen = (s.method_DB[keyword])[s._INDX_DATA_TO_GEN]
                        if (data_to_gen == []):
                            # In this case, FDPS directive for the clear method does not exist.
                            # In this case, we apply the default initialization to all the 
                            # member variables.
                            gen_mode = 0
                            # (1) Set members, data_types, array_dims, vals
                            members = []
                            data_types = []
                            array_dims = []
                            vals = []
                            for mbr in s.members:
                                fdpsDT, array_dim = self.__get_fdpsDT(mbr)
                                val = self.__default_IV_of_cppDT[fdpsDT]
                                members.append(mbr.name)
                                data_types.append(fdpsDT)
                                array_dims.append(array_dim)
                                vals.append(val)
                            # (2) Set method_DB
                            (s.method_DB[keyword])[s._INDX_GENERABLE_FLAG] = True
                            (s.method_DB[keyword])[s._INDX_DATA_TO_GEN] = [gen_mode,members,data_types,array_dims,vals]

        # Make the lists of FullParticle, EssentialParticleI, EssentialParticleJ, and Force
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    print("... {0}".format(s.name))
                    if (s.attrib["FP"] == True):
                        self.__FPs.append(s.name)
                    if (s.attrib["EPI"] == True):
                        self.__EPIs.append(s.name)
                    if (s.attrib["EPJ"] == True):
                        self.__EPJs.append(s.name)
                    if (s.attrib["Force"] == True):
                        self.__Forces.append(s.name)
        # Check if the lists above is empty or not
        if (len(self.__FPs) == 0):
            self.__print_error("No FullParticle exists!")
            sys.exit()
        if (len(self.__EPIs) == 0):
            self.__print_error("No EssentialParticleI exists!")
            sys.exit()
        if (len(self.__EPJs) == 0):
            self.__print_error("No EssentialParticleJ exists!")
            sys.exit()
        if (len(self.__Forces) == 0):
            self.__print_error("No Force exists!")
            sys.exit()
        
        # Compute self.__tree_kinds
        # Here, we must check the existence of the method getRSearch()
        # for some kinds of tree (FDPS specification)
        # (1) treeForForceLong
        candidates = list(itertools.product(["Long"], \
                                            self.__Forces, \
                                            self.__EPIs, \
                                            self.__EPJs, \
                                            self.__multipole_kinds))
        targets = ["MonopoleWithScatterSearch", \
                   "QuadrupoleWithScatterSearch", \
                   "MonopoleWithCutoff"]
        __tree_kinds_long = []
        for tree_t,force_t,epi_t,epj_t,multipole_t in candidates:
            item = [tree_t,force_t,epi_t,epj_t,multipole_t]
            if (multipole_t in targets):
                if (self.__getRSearch_existence_check(epj_t)):
                    __tree_kinds_long.append(item)
            else:
                __tree_kinds_long.append(item)
        # (2) treeForForceShort
        candidates = list(itertools.product(["Short"], \
                                             self.__Forces, \
                                             self.__EPIs, \
                                             self.__EPJs, \
                                             self.__neighbor_kinds))
        targets = ["Gather","Scatter","Symmetry"]
        __tree_kinds_short = []
        for tree_t,force_t,epi_t,epj_t,search_t in candidates:
            item = [tree_t,force_t,epi_t,epj_t,search_t]
            if (search_t in targets):
                if (self.__getRSearch_existence_check(epj_t)):
                    __tree_kinds_short.append(item)
            else:
                __tree_kinds_short.append(item)
        # (3) Pack
        self.__tree_kinds = __tree_kinds_long + __tree_kinds_short
        # [Option] check the content
        #print(len(self.__tree_kinds))
        #print(self.__tree_kinds)
        #sys.exit()

        # Compute self.__calc_force_kinds
        # (1) calcForceAll* for short-range interaction
        if (len(__tree_kinds_short) > 0):
            candidates = list(itertools.product(self.__FPs, \
                                                __tree_kinds_short))
            for fp_t,[tree_t,force_t,epi_t,epj_t,search_t] in candidates:
                item = [fp_t,tree_t,force_t,epi_t,epj_t,search_t]
                # (i) check if copyFromFP of epi_t copies from fp_t
                ret1 = self.__check_copyFromFP_consistency(epi_t,fp_t)
                # (ii) check if copyFromFP of epj_t copies from fp_t
                ret2 = self.__check_copyFromFP_consistency(epj_t,fp_t)
                # (iii) check if copyFromForce of fp_t copies from force_t
                ret3 = self.__check_copyFromForce_consistency(fp_t,force_t)
                if (ret1 and ret2 and ret3):
                    self.__calc_force_kinds_short.append(item)
        # (2) calcForceAll* for long-range interaction
        if (len(__tree_kinds_long) > 0):
            candidates = list(itertools.product(self.__FPs, \
                                                __tree_kinds_long))
            for fp_t,[tree_t,force_t,epi_t,epj_t,multipole_t] in candidates:
                item = [fp_t,tree_t,force_t,epi_t,epj_t,multipole_t]
                # (i) check if copyFromFP of epi_t copies from fp_t
                ret1 = self.__check_copyFromFP_consistency(epi_t,fp_t)
                # (ii) check if copyFromFP of epj_t copies from fp_t
                ret2 = self.__check_copyFromFP_consistency(epj_t,fp_t)
                # (iii) check if copyFromForce of fp_t copies from force_t
                ret3 = self.__check_copyFromForce_consistency(fp_t,force_t)
                if (ret1 and ret2 and ret3):
                    self.__calc_force_kinds_long.append(item)
        # [Option] check the content
        #print(len(self.__calc_force_kinds_short))
        #print(self.__calc_force_kinds_short)
        #print(len(self.__calc_force_kinds_long))
        #print(self.__calc_force_kinds_long)
        #sys.exit()

    def check(self):
        # Output to stdout
        self.__print_checkpoint("checking consistency....")

        #-----------------------------------------------------------------
        # [1] Check the consistency of method_DB
        for f in self.files:
            print("--> the file {0}".format(f.name))
            for mod in f.modules:
                print("-----> the module {0}".format(mod.name))
                for s in mod.structures:
                    print("--------> the structure {0}".format(s.name))
                    # Flag necessary methods
                    for key in s.attrib.keys():
                        if (s.attrib[key] == True):
                            if (key == "FP"):
                                s.method_DB["getPos"][s._INDX_NECESSARY_FLAG] = True
                                s.method_DB["copyFromForce"][s._INDX_NECESSARY_FLAG] = True
                            elif ((key == "EPI") or (key == "EPJ")):
                                s.method_DB["getPos"][s._INDX_NECESSARY_FLAG] = True
                                s.method_DB["copyFromFP"][s._INDX_NECESSARY_FLAG] = True
                            elif (key == "Force"):
                                s.method_DB["clear"][s._INDX_NECESSARY_FLAG] = True
                    # Check if all the necessary methods are generable
                    for key in s.method_DB.keys():
                        if (s.method_DB[key][s._INDX_NECESSARY_FLAG] == True): # if necessary method is True
                            if (s.method_DB[key][s._INDX_GENERABLE_FLAG] == False):
                                msg = "necessary method {0} cannot generate for structure {1}".format(key,s.name)
                                self.__print_error(msg)
                                sys.exit()
                    # Check if the class names in the arguments of copyFromForce and 
                    # copyFromFP really exist
                    keys = ["copyFromForce", "copyFromFP"]
                    for key in keys:
                        if (s.method_DB[key][s._INDX_GENERABLE_FLAG] == True):
                            meth_gen_data = (s.method_DB[key])[s._INDX_DATA_TO_GEN]
                            for data in meth_gen_data:
                                class_name = data.class_name
                                is_exist = self.__class_existence_check(class_name)
                                if (is_exist == False):
                                    self.__print_error("the structure {0} not found.".format(class_name))
                                    sys.exit()
                    # (here, we also check the consistency of member variables)
        
        # Check the result
        #for f in self.files:
        #    for mod in f.modules:
        #        # Output the module info.
        #        print("===============================================")
        #        print(" Module name  : {0}".format(mod.name))
        #        print(" buf_addr_sta : {0}".format(mod.buf_addr_sta))
        #        print(" buf_addr_end : {0}".format(mod.buf_addr_end))
        #        # Output the structure info.
        #        for s in mod.structures:
        #            print("-----------------------------------------------")
        #            print("< struct name : {0} >".format(s.name))
        #            for key,val in s.attrib.items():
        #                print("   {0} : {1}".format(key,val))
        #            print(" buf_addr_sta : {0}".format(s.buf_addr_sta))
        #            print(" buf_addr_end : {0}".format(s.buf_addr_end))
        #            print("** members")
        #            for mbr in s.members:
        #                if (mbr.is_array):
        #                    # In this case, this member is an array
        #                    print("{0}, dimension({1}) :: {2} {3}".format(mbr.data_type, \
        #                                                                  mbr.array_dim, \
        #                                                                  mbr.name, \
        #                                                                  mbr.attrib))
        #                else:
        #                    # In this case, this member is a scalar
        #                    print("{0} :: {1} {2}".format(mbr.data_type,mbr.name,mbr.attrib))

    def generate(self,output_dir):
        # Output to stdout
        self.__print_checkpoint("generating FDPS Fortran interface programs ....")
        # Make the directory
        output_dir = output_dir.strip()
        if (os.path.exists(output_dir) is False):
            os.mkdirs(output_dir)
        # Compute blueprint_dir 
        blueprint_dir = os.path.abspath(os.path.dirname(__file__)) \
                      + "/../src/fortran_interface/blueprints"
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
        # (2) user_defined.cpp
        file_name = output_dir + "/user_defined.hpp"
        fh = open(file_name,"w")
        # Write the header of the files
        fh.write("#pragma once\n")
        #fh.write("#include <cstdint>\n")
        fh.write("#include <cstring>\n") # To use std::strcpy()
        fh.write("#include <particle_simulator.hpp>\n\n")
        # Write the definitions of C++ classes
        for f in self.files:
            for mod in f.modules:
                for s in mod.structures:
                    # Write the header of the class
                    fh.write("class {0} {{\n".format(s.name))
                    fh.write("public:\n")
                    # Write the definitions of the member variables
                    for mbr in s.members:
                        # Get data type in FDPS-style
                        data_type, array_dim = self.__get_fdpsDT(mbr)
                        if (array_dim != 0):
                            fh.write("   {0} {1}[{2}];\n".format(data_type, \
                                                                 mbr.name,  \
                                                                 array_dim))
                        else:
                            fh.write("   {0} {1};\n".format(data_type, \
                                                            mbr.name))
                    fh.write("\n")
                    # Write the definitions of the method
                    if (s.method_DB["getId"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["getId"][s._INDX_DATA_TO_GEN]
                        self.__write_getId(fh,data_type,member_name)
                    if (s.method_DB["getPos"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["getPos"][s._INDX_DATA_TO_GEN]
                        self.__write_getPos(fh,data_type,member_name)
                    if (s.method_DB["setPos"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["setPos"][s._INDX_DATA_TO_GEN]
                        self.__write_setPos(fh,data_type,member_name)
                    if (s.method_DB["getCharge"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["getCharge"][s._INDX_DATA_TO_GEN]
                        self.__write_getCharge(fh,data_type,member_name)
                    if (s.method_DB["getChargePM"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["getChargePM"][s._INDX_DATA_TO_GEN]
                        self.__write_getChargePM(fh,data_type,member_name)
                    if (s.method_DB["getRSearch"][s._INDX_GENERABLE_FLAG] == True):
                        data_type, member_name = s.method_DB["getRSearch"][s._INDX_DATA_TO_GEN]
                        self.__write_getRSearch(fh,data_type,member_name)
                    if (s.method_DB["copyFromForce"][s._INDX_GENERABLE_FLAG] == True):
                        meth_gen_data = s.method_DB["copyFromForce"][s._INDX_DATA_TO_GEN]
                        self.__write_copyFromForce(fh,meth_gen_data)
                    if (s.method_DB["copyFromForcePM"][s._INDX_GENERABLE_FLAG] == True):
                        member_name = s.method_DB["copyFromForcePM"][s._INDX_DATA_TO_GEN]
                        self.__write_copyFromForcePM(fh,member_name)
                    if (s.method_DB["copyFromFP"][s._INDX_GENERABLE_FLAG] == True):
                        meth_gen_data = s.method_DB["copyFromFP"][s._INDX_DATA_TO_GEN]
                        self.__write_copyFromFP(fh,meth_gen_data)
                    if (s.method_DB["clear"][s._INDX_GENERABLE_FLAG] == True):
                        meth_gen_data = s.method_DB["clear"][s._INDX_DATA_TO_GEN]
                        self.__write_clear(fh,meth_gen_data)
                    # Write the footer of the class
                    fh.write("};\n\n")
                    # Write the implementation of the clear method if needed
                    if (s.method_DB["clear"][s._INDX_GENERABLE_FLAG] == True):
                        meth_gen_data = s.method_DB["clear"][s._INDX_DATA_TO_GEN]
                        force_name = s.name
                        self.__write_clear_impl(fh,meth_gen_data,force_name)
        fh.close()
        #---------------------------------------------------------------------------------
        # (3) FDPS_Manipulators.cpp
        # Pattern list
        pat_create_psys          = re.compile(r"fdps-autogen:create_psys;")
        pat_init_psys            = re.compile(r"fdps-autogen:init_psys;")
        pat_get_psys_memsize     = re.compile(r"fdps-autogen:get_psys_memsize;")
        pat_get_psys_time_prof   = re.compile(r"fdps-autogen:get_psys_time_prof;")
        pat_clear_psys_time_prof = re.compile(r"fdps-autogen:clear_psys_time_prof;")
        pat_set_nptcl_smpl       = re.compile(r"fdps-autogen:set_nptcl_smpl;")
        pat_set_nptcl_loc        = re.compile(r"fdps-autogen:set_nptcl_loc;")
        pat_get_nptcl_loc        = re.compile(r"fdps-autogen:get_nptcl_loc;")
        pat_get_nptcl_glb        = re.compile(r"fdps-autogen:get_nptcl_glb;")
        pat_get_psys_cptr        = re.compile(r"fdps-autogen:get_psys_cptr;")
        pat_exch_ptcl            = re.compile(r"fdps-autogen:exchange_particle;")
        pat_add_ptcl             = re.compile(r"fdps-autogen:add_particle;")
        pat_sort_ptcl            = re.compile(r"fdps-autogen:sort_particle;")
        pat_remove_ptcl          = re.compile(r"fdps-autogen:remove_particle;")
        pat_adjust_pos           = re.compile(r"fdps-autogen:adjust_pos_into_root_domain;")
        pat_col_smpl_ptcl        = re.compile(r"fdps-autogen:collect_sample_particle;")
        pat_dd_all               = re.compile(r"fdps-autogen:decompose_domain_all;")
        pat_create_tree          = re.compile(r"fdps-autogen:create_tree;")
        pat_init_tree            = re.compile(r"fdps-autogen:init_tree;")
        pat_get_tree_memsize     = re.compile(r"fdps-autogen:get_tree_memsize;")
        pat_get_tree_time_prof   = re.compile(r"fdps-autogen:get_tree_time_prof;")
        pat_clear_tree_time_prof = re.compile(r"fdps-autogen:clear_tree_time_prof;")
        pat_get_nint_ep_ep_loc   = re.compile(r"fdps-autogen:get_num_interact_ep_ep_loc;")
        pat_get_nint_ep_sp_loc   = re.compile(r"fdps-autogen:get_num_interact_ep_sp_loc;")
        pat_get_nint_ep_ep_glb   = re.compile(r"fdps-autogen:get_num_interact_ep_ep_glb;")
        pat_get_nint_ep_sp_glb   = re.compile(r"fdps-autogen:get_num_interact_ep_sp_glb;")
        pat_clear_nint           = re.compile(r"fdps-autogen:clear_num_interact;")
        pat_get_nwalk_loc        = re.compile(r"fdps-autogen:get_num_tree_walk_loc;")
        pat_get_nwalk_glb        = re.compile(r"fdps-autogen:get_num_tree_walk_glb;")
        pat_calc_force_all_WB    = re.compile(r"fdps-autogen:calc_force_all_and_write_back;")
        pat_calc_force_all       = re.compile(r"fdps-autogen:calc_force_all;")
        pat_calc_force_MT        = re.compile(r"fdps-autogen:calc_force_making_tree;")
        pat_calc_force_WB        = re.compile(r"fdps-autogen:calc_force_and_write_back;")
        pat_get_ngb_list         = re.compile(r"fdps-autogen:get_neighbor_list;")
        pat_get_epj_from_id      = re.compile(r"fdps-autogen:get_epj_from_id;")
        pat_set_psys_pm          = re.compile(r"fdps-autogen:set_psys_of_pm;")
        pat_pm_force_all_WB      = re.compile(r"fdps-autogen:calc_pm_force_all_and_write_back;")
        # Auto-generation
        input_file  = blueprint_dir + "/FDPS_Manipulators_blueprint.cpp"
        output_file = output_dir    + "/FDPS_Manipulators.cpp"
        ifh = open(input_file,'r')
        ofh = open(output_file,'w')
        for line in ifh:
            ofh.write(line) # simply copy
            #----- psys APIs ------
            if (pat_create_psys.search(line)):
                self.__write_psys_branch_for_create_psys(ofh)
            if (pat_init_psys.search(line)):
                inst = "psys->initialize();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_get_psys_memsize.search(line)):
                inst = "return (long long int) psys->getMemSizeUsed();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_get_psys_time_prof.search(line)):
                inst = "*prof = psys->getTimeProfile();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_clear_psys_time_prof.search(line)):
                inst = "psys->clearTimeProfile();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_set_nptcl_smpl.search(line)):
                inst = "psys->setAverageTargetNumberOfSampleParticlePerProcess(numPtcl);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_set_nptcl_loc.search(line)):
                inst = "psys->setNumberOfParticleLocal(numPtcl);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nptcl_loc.search(line)):
                inst = "return psys->getNumberOfParticleLocal();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nptcl_glb.search(line)):
                inst = "return psys->getNumberOfParticleGlobal();"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_get_psys_cptr.search(line)):
                inst = "*cptr = (void *) &((*psys)[0]);"
                cmnt = "// [!!IMPORTANT!!] " \
                     + "The operator [0] is necessary to obtain the correct address of " \
                     + "the particle system object." 
                self.__write_psys_branch_for_normal_APIs(ofh,inst,cmnt)
            if (pat_exch_ptcl.search(line)):
                inst = "psys->exchangeParticle(*dinfo);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_add_ptcl.search(line)):
                self.__write_add_particle_cpp_impl(ofh)
            if (pat_sort_ptcl.search(line)):
                self.__write_sort_particle_cpp_impl(ofh)
            if (pat_remove_ptcl.search(line)):
                inst = "psys->removeParticle(ptcl_indx,numPtcl);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_adjust_pos.search(line)):
                inst = "psys->adjustPositionIntoRootDomain(*dinfo);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            #----- dinfo APIs ------
            if (pat_col_smpl_ptcl.search(line)):
                inst = "dinfo->collectSampleParticle(*psys,clear,weight);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_dd_all.search(line)):
                inst = "dinfo->decomposeDomainAll(*psys,weight);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            #----- tree APIs ------
            if (pat_create_tree.search(line)):
                self.__write_tree_branch_for_create_tree(ofh)
            if (pat_init_tree.search(line)):
                inst = "tree->initialize(numPtcl,theta,n_leaf_limit,n_group_limit);"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_tree_memsize.search(line)):
                inst = "return (long long int) tree->getMemSizeUsed();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_tree_time_prof.search(line)):
                inst = "*prof = tree->getTimeProfile();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_clear_tree_time_prof.search(line)):
                inst = "tree->clearTimeProfile();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nint_ep_ep_loc.search(line)):
                inst = "return (long long int) tree->getNumberOfInteractionEPEPLocal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nint_ep_sp_loc.search(line)):
                inst = "return (long long int) tree->getNumberOfInteractionEPSPLocal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nint_ep_ep_glb.search(line)):
                inst = "return (long long int) tree->getNumberOfInteractionEPEPGlobal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nint_ep_sp_glb.search(line)):
                inst = "return (long long int) tree->getNumberOfInteractionEPSPGlobal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_clear_nint.search(line)):
                inst = "tree->clearNumberOfInteraction();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nwalk_loc.search(line)):
                inst = "return (long long int) tree->getNumberOfWalkLocal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            if (pat_get_nwalk_glb.search(line)):
                inst = "return (long long int) tree->getNumberOfWalkGlobal();"
                self.__write_tree_branch_for_normal_APIs(ofh,inst)
            #----- tree APIs (calc_force*) ------
            if (pat_calc_force_all_WB.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                fdps_method_name  = "calcForceAllAndWriteBack"
                self.__write_calc_force_branch(ofh,generic_func_name,fdps_method_name)
            if (pat_calc_force_all.search(line)):
                generic_func_name = "calc_force_all"
                fdps_method_name  = "calcForceAll"
                self.__write_calc_force_branch(ofh,generic_func_name,fdps_method_name)
            if (pat_calc_force_MT.search(line)):
                generic_func_name = "calc_force_making_tree"
                fdps_method_name  = "calcForceMakingTree"
                self.__write_calc_force_branch(ofh,generic_func_name,fdps_method_name)
            if (pat_calc_force_WB.search(line)):
                generic_func_name = "calc_force_and_write_back"
                fdps_method_name  = "calcForceAndWriteBack"
                self.__write_calc_force_branch(ofh,generic_func_name,fdps_method_name)
            #----- tree APIs (get_neighbor_list & get_epj_from_id*) ------
            if (pat_get_ngb_list.search(line)):
                self.__write_get_ngb_list_branch(ofh)
            if (pat_get_epj_from_id.search(line)):
                self.__write_get_epj_from_id_branch(ofh)
            #----- ParticleMesh APIs ------
            if (pat_set_psys_pm.search(line)):
                inst = "pm->setParticleParticleMesh(*psys,clear);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
            if (pat_pm_force_all_WB.search(line)):
                inst = "pm->calcForceAllAndWriteBack(*psys,*dinfo);"
                self.__write_psys_branch_for_normal_APIs(ofh,inst)
        ifh.close()
        ofh.close()
        #---------------------------------------------------------------------------------
        # (4) FDPS_Manipulators.h 
        # Pattern list
        pat_add_ptcl          = re.compile(r"fdps-autogen:add_particle;")
        pat_sort_ptcl         = re.compile(r"fdps-autogen:sort_particle;")
        pat_calc_force_all_WB = re.compile(r"fdps-autogen:calc_force_all_and_write_back;")
        pat_calc_force_all    = re.compile(r"fdps-autogen:calc_force_all;")
        pat_calc_force_MT     = re.compile(r"fdps-autogen:calc_force_making_tree;")
        pat_calc_force_WB     = re.compile(r"fdps-autogen:calc_force_and_write_back;")
        pat_get_ngb_list      = re.compile(r"fdps-autogen:get_neighbor_list;")
        pat_get_epj_from_id   = re.compile(r"fdps-autogen:get_epj_from_id;")
        # Auto-generation
        input_file  = blueprint_dir + "/FDPS_Manipulators_blueprint.h"
        output_file = output_dir    + "/FDPS_Manipulators.h"
        ifh = open(input_file,'r')
        ofh = open(output_file,'w')
        for line in ifh:
            ofh.write(line) # simply copy
            #----- psys APIs -----
            if (pat_add_ptcl.search(line)):
                self.__write_add_particle_protoDCL(ofh)
            if (pat_sort_ptcl.search(line)):
                self.__write_sort_particle_protoDCL(ofh)
            #----- tree APIs (calc_force*) -----
            if (pat_calc_force_all_WB.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                self.__write_calc_force_protoDCL_branch(ofh,generic_func_name)
            if (pat_calc_force_all.search(line)):
                generic_func_name = "calc_force_all"
                self.__write_calc_force_protoDCL_branch(ofh,generic_func_name)
            if (pat_calc_force_MT.search(line)):
                generic_func_name = "calc_force_making_tree"
                self.__write_calc_force_protoDCL_branch(ofh,generic_func_name)
            if (pat_calc_force_WB.search(line)):
                generic_func_name = "calc_force_and_write_back"
                self.__write_calc_force_protoDCL_branch(ofh,generic_func_name)
            #----- tre APIs (get_neighbor_list* & get_epj_from_id*)
            if (pat_get_ngb_list.search(line)):
                self.__write_get_ngb_list_protoDCL_branch(ofh)
            if (pat_get_epj_from_id.search(line)):
                self.__write_get_epj_from_id_protoDCL_branch(ofh)
        ifh.close()
        ofh.close()
        #---------------------------------------------------------------------------------
        # (5) FDPS_ftn_if.cpp
        # Pattern list
        pat_add_ptcl          = re.compile(r"fdps-autogen:add_particle")
        pat_sort_ptcl         = re.compile(r"fdps-autogen:sort_particle")
        pat_calc_force_all_WB = re.compile(r"fdps-autogen:calc_force_all_and_write_back;")
        pat_calc_force_all    = re.compile(r"fdps-autogen:calc_force_all;")
        pat_calc_force_MT     = re.compile(r"fdps-autogen:calc_force_making_tree;")
        pat_calc_force_WB     = re.compile(r"fdps-autogen:calc_force_and_write_back;")
        pat_get_ngb_list      = re.compile(r"fdps-autogen:get_neighbor_list;")
        pat_get_epj_from_id   = re.compile(r"fdps-autogen:get_epj_from_id;")
        pat_get_min_value     = re.compile(r"fdps-autogen:get_min_value;")
        pat_get_max_value     = re.compile(r"fdps-autogen:get_max_value;")
        pat_get_sum           = re.compile(r"fdps-autogen:get_sum;")
        pat_broadcast         = re.compile(r"fdps-autogen:broadcast;")
        # Auto-generation
        input_file  = blueprint_dir + "/FDPS_ftn_if_blueprint.cpp"
        output_file = output_dir    + "/FDPS_ftn_if.cpp"
        ifh = open(input_file,'r')
        ofh = open(output_file,'w')
        for line in ifh:
            ofh.write(line) # simply copy
            #----- psys APIs -----
            if (pat_add_ptcl.search(line)):
                self.__write_add_particle_cpp_if(ofh)
            if (pat_sort_ptcl.search(line)):
                self.__write_sort_particle_cpp_if(ofh)
            #----- tree APIs -----
            if (pat_calc_force_all_WB.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                self.__write_calc_force_cpp_if(ofh,generic_func_name)
            if (pat_calc_force_all.search(line)):
                generic_func_name = "calc_force_all"
                self.__write_calc_force_cpp_if(ofh,generic_func_name)
            if (pat_calc_force_MT.search(line)):
                generic_func_name = "calc_force_making_tree"
                self.__write_calc_force_cpp_if(ofh,generic_func_name)
            if (pat_calc_force_WB.search(line)):
                generic_func_name = "calc_force_and_write_back"
                self.__write_calc_force_cpp_if(ofh,generic_func_name)
            if (pat_get_ngb_list.search(line)):
                self.__write_get_ngb_list_cpp_if(ofh)
            if (pat_get_epj_from_id.search(line)):
                self.__write_get_epj_from_id_cpp_if(ofh)
            #----- comm APIs -----
            if (pat_get_min_value.search(line)):
                generic_func_name = "get_min_value"
                fdps_api_name = "PS::Comm::getMinValue"
                self.__write_reduction_cpp_if(ofh,generic_func_name,fdps_api_name,True)
            if (pat_get_max_value.search(line)):
                generic_func_name = "get_max_value"
                fdps_api_name = "PS::Comm::getMaxValue"
                self.__write_reduction_cpp_if(ofh,generic_func_name,fdps_api_name,True)
            if (pat_get_sum.search(line)):
                generic_func_name = "get_sum"
                fdps_api_name = "PS::Comm::getSum"
                self.__write_reduction_cpp_if(ofh,generic_func_name,fdps_api_name,False)
            if (pat_broadcast.search(line)):
                self.__write_broadcast_cpp_if(ofh)
        ofh.close()
        ifh.close()
        #---------------------------------------------------------------------------------
        # (6) FDPS_module.F90
        # Pattern list
        pat_get_psys_fptr_meth       = re.compile(r"fdps-autogen:get_psys_fptr:method;")
        pat_get_psys_fptr_decl       = re.compile(r"fdps-autogen:get_psys_fptr:decl;")
        pat_get_psys_fptr_impl       = re.compile(r"fdps-autogen:get_psys_fptr:impl;")
        pat_add_ptcl_meth            = re.compile(r"fdps-autogen:add_particle:method;")
        pat_add_ptcl_decl            = re.compile(r"fdps-autogen:add_particle:decl;")
        pat_add_ptcl_if              = re.compile(r"fdps-autogen:add_particle:if")
        pat_add_ptcl_impl            = re.compile(r"fdps-autogen:add_particle:impl")
        pat_sort_ptcl_if             = re.compile(r"fdps-autogen:sort_particle:if")
        pat_sort_ptcl_impl           = re.compile(r"fdps-autogen:sort_particle:impl")
        pat_calc_force_all_WB_if     = re.compile(r"fdps-autogen:calc_force_all_and_write_back:if;")
        pat_calc_force_all_WB_impl_s = re.compile(r"fdps-autogen:calc_force_all_and_write_back:impl:short;") 
        pat_calc_force_all_WB_impl_l = re.compile(r"fdps-autogen:calc_force_all_and_write_back:impl:long;")
        pat_calc_force_all_if        = re.compile(r"fdps-autogen:calc_force_all:if;")
        pat_calc_force_all_impl_s    = re.compile(r"fdps-autogen:calc_force_all:impl:short;") 
        pat_calc_force_all_impl_l    = re.compile(r"fdps-autogen:calc_force_all:impl:long;")
        pat_calc_force_MT_if         = re.compile(r"fdps-autogen:calc_force_making_tree:if;")
        pat_calc_force_MT_impl_s     = re.compile(r"fdps-autogen:calc_force_making_tree:impl:short;") 
        pat_calc_force_MT_impl_l     = re.compile(r"fdps-autogen:calc_force_making_tree:impl:long;")
        pat_calc_force_WB_if         = re.compile(r"fdps-autogen:calc_force_and_write_back:if;")
        pat_calc_force_WB_impl_s     = re.compile(r"fdps-autogen:calc_force_and_write_back:impl:short;") 
        pat_calc_force_WB_impl_l     = re.compile(r"fdps-autogen:calc_force_and_write_back:impl:long;")
        pat_get_ngb_list_meth        = re.compile(r"fdps-autogen:get_neighbor_list:method;")
        pat_get_ngb_list_decl        = re.compile(r"fdps-autogen:get_neighbor_list:decl;")
        pat_get_ngb_list_if          = re.compile(r"fdps-autogen:get_neighbor_list:if;")
        pat_get_ngb_list_impl        = re.compile(r"fdps-autogen:get_neighbor_list:impl;")
        pat_get_epj_from_id_meth     = re.compile(r"fdps-autogen:get_epj_from_id:method;")
        pat_get_epj_from_id_decl     = re.compile(r"fdps-autogen:get_epj_from_id:decl;")
        pat_get_epj_from_id_if       = re.compile(r"fdps-autogen:get_epj_from_id:if;")
        pat_get_epj_from_id_impl     = re.compile(r"fdps-autogen:get_epj_from_id:impl;")
        pat_get_min_value_meth       = re.compile(r"fdps-autogen:get_min_value:method;")
        pat_get_min_value_decl       = re.compile(r"fdps-autogen:get_min_value:decl;")
        pat_get_min_value_ftn_if     = re.compile(r"fdps-autogen:get_min_value:ftn_if;")
        pat_get_min_value_ftn_impl   = re.compile(r"fdps-autogen:get_min_value:ftn_impl;")
        pat_get_max_value_meth       = re.compile(r"fdps-autogen:get_max_value:method;")
        pat_get_max_value_decl       = re.compile(r"fdps-autogen:get_max_value:decl;")
        pat_get_max_value_ftn_if     = re.compile(r"fdps-autogen:get_max_value:ftn_if;")
        pat_get_max_value_ftn_impl   = re.compile(r"fdps-autogen:get_max_value:ftn_impl;")
        pat_get_sum_meth             = re.compile(r"fdps-autogen:get_sum:method;")
        pat_get_sum_decl             = re.compile(r"fdps-autogen:get_sum:decl;")
        pat_get_sum_ftn_if           = re.compile(r"fdps-autogen:get_sum:ftn_if;")
        pat_get_sum_ftn_impl         = re.compile(r"fdps-autogen:get_sum:ftn_impl;")
        pat_broadcast_meth           = re.compile(r"fdps-autogen:broadcast:method;")
        pat_broadcast_decl           = re.compile(r"fdps-autogen:broadcast:decl;")
        pat_broadcast_ftn_if         = re.compile(r"fdps-autogen:broadcast:ftn_if;")
        pat_broadcast_ftn_impl       = re.compile(r"fdps-autogen:broadcast:ftn_impl;")
        # Auto-generation
        input_file  = blueprint_dir + "/FDPS_module_blueprint.F90"
        output_file = output_dir    + "/FDPS_module.F90"
        ifh = open(input_file,'r')
        ofh = open(output_file,'w')
        for line in ifh:
            ofh.write(line) # simply copy
            #----- psys APIs -----
            if (pat_get_psys_fptr_meth.search(line)):
                self.__write_get_psys_fptr_meth(ofh)
            if (pat_get_psys_fptr_decl.search(line)):
                self.__write_get_psys_fptr_decl(ofh)
            if (pat_get_psys_fptr_impl.search(line)):
                self.__write_get_psys_fptr_impl(ofh)

            if (pat_add_ptcl_meth.search(line)):
                self.__write_add_particle_meth(ofh)
            if (pat_add_ptcl_decl.search(line)):
                self.__write_add_particle_decl(ofh)
            if (pat_add_ptcl_if.search(line)):
                self.__write_add_particle_ftn_if(ofh)
            if (pat_add_ptcl_impl.search(line)):
                self.__write_add_particle_impl(ofh)

            if (pat_sort_ptcl_if.search(line)):
                self.__write_sort_particle_ftn_if(ofh)
            if (pat_sort_ptcl_impl.search(line)):
                self.__write_sort_particle_impl(ofh)
            #----- tree APIs (calc_force) -----
            # (1) calc_force_all_and_write_back
            if (pat_calc_force_all_WB_if.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                self.__write_calc_force_ftn_if(ofh,generic_func_name)
            if (pat_calc_force_all_WB_impl_s.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                self.__write_calc_force_ftn_impl_s(ofh,generic_func_name)
            if (pat_calc_force_all_WB_impl_l.search(line)):
                generic_func_name = "calc_force_all_and_write_back"
                self.__write_calc_force_ftn_impl_l(ofh,generic_func_name)
            # (2) calc_force_all
            if (pat_calc_force_all_if.search(line)):
                generic_func_name = "calc_force_all"
                self.__write_calc_force_ftn_if(ofh,generic_func_name)
            if (pat_calc_force_all_impl_s.search(line)):
                generic_func_name = "calc_force_all"
                self.__write_calc_force_ftn_impl_s(ofh,generic_func_name)
            if (pat_calc_force_all_impl_l.search(line)):
                generic_func_name = "calc_force_all"
                self.__write_calc_force_ftn_impl_l(ofh,generic_func_name)
            # (3) calc_force_making_tree
            if (pat_calc_force_MT_if.search(line)):
                generic_func_name = "calc_force_making_tree"
                self.__write_calc_force_ftn_if(ofh,generic_func_name)
            if (pat_calc_force_MT_impl_s.search(line)):
                generic_func_name = "calc_force_making_tree"
                self.__write_calc_force_ftn_impl_s(ofh,generic_func_name)
            if (pat_calc_force_MT_impl_l.search(line)):
                generic_func_name = "calc_force_making_tree"
                self.__write_calc_force_ftn_impl_l(ofh,generic_func_name)
            # (4) calc_force_and_write_back
            if (pat_calc_force_WB_if.search(line)):
                generic_func_name = "calc_force_and_write_back"
                self.__write_calc_force_ftn_if(ofh,generic_func_name)
            if (pat_calc_force_WB_impl_s.search(line)):
                generic_func_name = "calc_force_and_write_back"
                self.__write_calc_force_ftn_impl_s(ofh,generic_func_name)
            if (pat_calc_force_WB_impl_l.search(line)):
                generic_func_name = "calc_force_and_write_back"
                self.__write_calc_force_ftn_impl_l(ofh,generic_func_name)
            #----- tree APIs (get_neighbor_list* & get_epj_from_id*) ------
            # (1) get_neighbor_list
            if (pat_get_ngb_list_meth.search(line)):
                self.__write_get_ngb_list_meth(ofh)
            if (pat_get_ngb_list_decl.search(line)):
                self.__write_get_ngb_list_decl(ofh)
            if (pat_get_ngb_list_if.search(line)):
                self.__write_get_ngb_list_ftn_if(ofh)
            if (pat_get_ngb_list_impl.search(line)):
                self.__write_get_ngb_list_impl(ofh)
            # (2) get_epj_from_id
            if (pat_get_epj_from_id_meth.search(line)):
                self.__write_get_epj_from_id_meth(ofh)
            if (pat_get_epj_from_id_decl.search(line)):
                self.__write_get_epj_from_id_decl(ofh)
            if (pat_get_epj_from_id_if.search(line)):
                self.__write_get_epj_from_id_ftn_if(ofh)
            if (pat_get_epj_from_id_impl.search(line)):
                self.__write_get_epj_from_id_impl(ofh)
            #----- comm APIs -----
            # (1) get_min_value
            if (pat_get_min_value_meth.search(line)):
                generic_func_name = "get_min_value"
                self.__write_reduction_meth(ofh,generic_func_name,True)
            if (pat_get_min_value_decl.search(line)):
                generic_func_name = "get_min_value"
                self.__write_reduction_decl(ofh,generic_func_name,True)
            if (pat_get_min_value_ftn_if.search(line)):
                generic_func_name = "get_min_value"
                self.__write_reduction_ftn_if(ofh,generic_func_name,True)
            if (pat_get_min_value_ftn_impl.search(line)):
                generic_func_name = "get_min_value"
                self.__write_reduction_ftn_impl(ofh,generic_func_name,True)
            # (2) get_max_value
            if (pat_get_max_value_meth.search(line)):
                generic_func_name = "get_max_value"
                self.__write_reduction_meth(ofh,generic_func_name,True)
            if (pat_get_max_value_decl.search(line)):
                generic_func_name = "get_max_value"
                self.__write_reduction_decl(ofh,generic_func_name,True)
            if (pat_get_max_value_ftn_if.search(line)):
                generic_func_name = "get_max_value"
                self.__write_reduction_ftn_if(ofh,generic_func_name,True)
            if (pat_get_max_value_ftn_impl.search(line)):
                generic_func_name = "get_max_value"
                self.__write_reduction_ftn_impl(ofh,generic_func_name,True)
            # (3) get_sum
            if (pat_get_sum_meth.search(line)):
                generic_func_name = "get_sum"
                self.__write_reduction_meth(ofh,generic_func_name,False)
            if (pat_get_sum_decl.search(line)):
                generic_func_name = "get_sum"
                self.__write_reduction_decl(ofh,generic_func_name,False)
            if (pat_get_sum_ftn_if.search(line)):
                generic_func_name = "get_sum"
                self.__write_reduction_ftn_if(ofh,generic_func_name,False)
            if (pat_get_sum_ftn_impl.search(line)):
                generic_func_name = "get_sum"
                self.__write_reduction_ftn_impl(ofh,generic_func_name,False)
            # (4) broadcast
            if (pat_broadcast_meth.search(line)):
                self.__write_broadcast_meth(ofh)
            if (pat_broadcast_decl.search(line)):
                self.__write_broadcast_decl(ofh)
            if (pat_broadcast_ftn_if.search(line)):
                self.__write_broadcast_ftn_if(ofh)
            if (pat_broadcast_ftn_impl.search(line)):
                self.__write_broadcast_ftn_impl(ofh)
        ofh.close()
        ifh.close()
        sys.exit()
        #---------------------------------------------------------------------------------

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
    description = "Analyze user's Fortran codes and generate C++/Fortran source files " \
                + "required to use FDPS from the user's Fortran code."
    parser = argparse.ArgumentParser(description=description)
    # Input option
    help_msg = "The PATHs of input Fortran files"
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

