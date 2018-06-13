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
    print("Module os,sys,shutil,struct,math are not found.")
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
class Automatic_Tester:
    def __init__(self):
        # These private variables stores the numbers of OpenMP threads
        # and MPI processes that are used for tests.
        self.__NTHRDS = 2
        self.__NPROCS = 2

        # The sample code `water` can be run with several modes.
        # self.__mode_info stores the information about each mode.
        # The content MUST BE consistent with Makefile at the top directory of `water`.
        self.__mode_info = collections.OrderedDict()
        self.__mode_info["rigid"] = {"src_dir":"./src/rigid",
                                     "mk_tgt_name":"",
                                     "exec_file_name":"main.out"}
        self.__mode_info["flex"]  = {"src_dir":"./src/flex",
                                     "mk_tgt_name":"",
                                     "exec_file_name":"main.out"}
        self.__mode_info["gpu"]   = {"src_dir":"./src/cuda",
                                     "mk_tgt_name":"gpu",
                                     "exec_file_name":"gpu.out"}
        self.__mode_info["reuse"] = {"src_dir":"./src/cuda",
                                     "mk_tgt_name":"reuse",
                                     "exec_file_name":"reuse.out"}

        # This private variable stores the original Makefile.
        self.__mkfile_org = collections.OrderedDict()
        self.__mkfile_org["rigid"] = []
        self.__mkfile_org["flex"]  = []
        self.__mkfile_org["gpu"]   = []
        self.__mkfile_org["reuse"] = []

        # This private variable stores some required files
        self.__depfile = collections.OrderedDict()
        self.__depfile["ice_unit.pdb"] = []

        # This private variable stores the values of CXX, CXXFLAGS,
        # in the original Makefile.
        self.__mkvars = []

        # This private variable stores an acceptable energy error
        self.__engy_err_crit = 1.0e-2

    def create_tests(self):
        # Set test configurations
        candidates = []
        # (1) w/o MPI
        CXX = ["time g++"]
        CXXFLAGS = ["", \
                    "-DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp"]
        RUNOPT = ["export OMP_NUM_THREADS={0}; ".format(self.__NTHRDS)]
        listtmp = self.__get_dictlist(CXX,CXXFLAGS,RUNOPT)
        candidates.extend(listtmp)
        # (2) w/ MPI
        CXX = ["time mpicxx"]
        CXXFLAGS = ["-DPARTICLE_SIMULATOR_MPI_PARALLEL", \
                    "-DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp"]
        RUNOPT = ["export OMP_NUM_THREADS={0}; mpirun --mca btl ^openib -np {1} ".format(self.__NTHRDS,self.__NPROCS)]
        listtmp = self.__get_dictlist(CXX,CXXFLAGS,RUNOPT)
        candidates.extend(listtmp)
        # Remove exceptions
        for i in range(0,len(candidates)):
            item = candidates[i]
            self.__mkvars.append(item) 
        # Check [for debug]
        #for item in self.__mkvars:
        #    print(item)

    def read_mkfiles(self):
        # Read Makefile in each mode
        for key,val in self.__mode_info.items():
            # Set filename
            fname = val["src_dir"] + "/Makefile"
            # Open file
            fp = open(fname,'r')
            # Read the file
            for line in fp:
                # Delete LF
                #line = line.rstrip()
                # Save the line
                self.__mkfile_org[key].append(line)
            # Close the file
            fp.close()
            # Check [for debug]
            #print(self.__mkfile_org[key])

    def read_depfiles(self):
        # Read required files
        for key,val in self.__depfile.items():
           # Set filename
           fname = key
           # Open file
           fp = open(fname,'r')
           # Read the file
           f = []
           for line in fp:
               # Delete LF
               #line = line.rstrip()
               # Save the line
               f.append(line)
           # Close the file
           fp.close()
           # Set
           self.__depfile[key] = f
           # Check [for debug]
           #print(self.__depfile[key])

    def __write_depfiles(self):
        for fname in self.__depfile.keys():
            fp = open(fname,'w')
            for line in self.__depfile[fname]:
                fp.write(line)
            fp.close()

    def __delete_depfiles(self):
        for fname in self.__depfile.keys():
            if os.path.isfile(fname):
                os.remove(fname)

    def test(self):
        # Get the current directory path
        root_dir = os.getcwd()
        print("root_dir = {0}".format(root_dir))
        # Test each mode
        for key,val in self.__mode_info.items():
            # Define various constants and declare lists
            pat_set_vars = re.compile(r"^\s*#\s*fdps-autotest-set-vars")
            pat_run_inst = re.compile(r"^\s*#\s*fdps-autotest-run")
            mkfile_for_test  = "Makefile.test"
            logfile_for_test = "stdout-test.txt"
            results = []
            # Perform test
            for mkvars in self.__mkvars:
                # Change the working directory
                os.chdir(root_dir) # just in case
                os.chdir(val["src_dir"])
                # Prepare required files
                self.__write_depfiles()
                # Delete Makefile.test if it exists
                if os.path.isfile(mkfile_for_test):
                    os.remove(mkfile_for_test)
                # Create a new Makefile
                fp = open(mkfile_for_test,'w')
                # Write to the file
                for line in self.__mkfile_org[key]:
                    # Copy
                    fp.write(line)
                    # Overwrite the makefile variables if this line
                    # corresponds to the insertion point.
                    if (pat_set_vars.match(line)):
                        text = """
                        #---------------------------------------------------
                        # !!! Overwrite the makefile variable for test !!!
                        CXX = {CXX}
                        CXXFLAGS = $(CXXFLAGS_COMMON) {CXXFLAGS}
                        #---------------------------------------------------
                        """.format(CXX=mkvars["CXX"], \
                                   CXXFLAGS=mkvars["CXXFLAGS"])
                        text = text[1:].rstrip()
                        text = textwrap.dedent(text)
                        text = text + "\n"
                        fp.write(text)
                    # Describe how to perform the test if this line
                    # corresponds to the insertion point.
                    if (pat_run_inst.match(line)):
                        text = """
                        #---------------------------------------------------
                        run:
                        \t{RUNOPT} ./{EXEC_FILE} > {LOGFILE} 2>&1
                        #---------------------------------------------------
                        """.format(RUNOPT=mkvars["RUNOPT"], \
                                   EXEC_FILE=val["exec_file_name"], \
                                   LOGFILE=logfile_for_test)
                        text = text[1:].rstrip()
                        text = textwrap.dedent(text)
                        text = text + "\n"
                        fp.write(text)
                # Close the file
                fp.close()
                # Make
                cmd = "make {0} -B -f {1}; echo $?".format(val["mk_tgt_name"],mkfile_for_test)
                outs,errs = subprocess.Popen(cmd, \
                                             stdout=subprocess.PIPE, \
                                             stderr=subprocess.PIPE, \
                                             shell=True).communicate()
                if (outs.rstrip() == "2"): # see https://linuxjm.osdn.jp/html/GNU_make/man1/make.1.html
                    msg = self.__get_red_text("failed to make. Stop this test.")
                    print("{0}".format(msg))
                    print("{0}".format(errs))
                    cmd = "make distclean -f {0}".format(mkfile_for_test)
                    self.__end_processing(results,cmd,[mkfile_for_test])
                    continue # Move to next test
                # Run
                cmd = "make run -f {0}".format(mkfile_for_test)
                subprocess.call(cmd,shell=True)
                # Obtain a string representing relative energy error
                cmd = "grep -e \"Rel.Err.\" {0} | awk 'END {{print $NF}}'".format(logfile_for_test)
                engy_err = subprocess.Popen(cmd, \
                                            stdout=subprocess.PIPE, \
                                            shell=True).communicate()[0]
                # Delete all kinds of blank characters from engy_err
                engy_err = engy_err.rstrip()    
                engy_err = engy_err.strip()     
                engy_err = engy_err.strip("\t") 
                # Error handling for engy_err
                if engy_err == "":
                    msg = self.__get_red_text("the program does not stop successfully.")
                    print("{0}".format(msg))
                    cmd = "make distclean -f {0}".format(mkfile_for_test)
                    self.__end_processing(results,cmd,[mkfile_for_test,logfile_for_test])
                    continue # Move to next test
                # Convert a string to a FP value 
                try:
                    engy_err = abs(float(engy_err))
                    #print(engy_err)
                except ValueError:
                    msg = self.__get_red_text("engy_err is not a FP value.")
                    print("{0}".format(msg))
                    cmd = "make distclean -f {0}".format(mkfile_for_test)
                    self.__end_processing(results,cmd,[mkfile_for_test,logfile_for_test])
                    continue # Move to next test
                # Store the result
                if (engy_err < self.__engy_err_crit):
                    results.append(0)
                else:
                    results.append(1)
                # Clean
                cmd = "make distclean -f {0}".format(mkfile_for_test)
                subprocess.call(cmd,shell=True)
                # Delete unnecessary files
                self.__delete_depfiles()
                os.remove(mkfile_for_test)
                os.remove(logfile_for_test)
            # Output the summary of the tests
            text = """
            #################################################################
               Summary of tests for sample/c++/water (run mode: {MODE})
            #################################################################
            """.format(MODE=key)
            text = text[1:].rstrip()
            text = textwrap.dedent(text)
            print(text)
            for testnum in range(0,len(self.__mkvars)):
                mkvars = self.__mkvars[testnum]
                result = results[testnum]
                if result == 0:
                    msg = self.__get_green_text("OK")
                else:
                    msg = self.__get_red_text("FAIL")
                text = """
                {TESTNUM}) CXX={CXX}, CXXFLAGS={CXXFLAGS} : {STATUS}
                """.format(TESTNUM=testnum, \
                           CXX=mkvars["CXX"], \
                           CXXFLAGS=mkvars["CXXFLAGS"], \
                           STATUS=msg)
                text = text[1:].rstrip()
                text = textwrap.dedent(text)
                print("{0}".format(text))
            # [DEBUG]
            if (key == "reuse"):
                sys.exit(0)
        # Check if all the tests are passed or not.
        return sum(results)

    def __get_dictlist(self,CXX,CXXFLAGS,RUNOPT):
        listtmp = list(itertools.product(CXX, \
                                         CXXFLAGS, \
                                         RUNOPT))
        dictlist = []
        for cxx,cxxflags,runopt in listtmp:
            item = collections.OrderedDict()
            item["CXX"] = cxx
            item["CXXFLAGS"] = cxxflags
            item["RUNOPT"] = runopt
            dictlist.append(item)
        return dictlist

    def __end_processing(self,results, \
                         cmd_to_be_run, \
                         files_to_be_removed):
        results.append(1)
        # Command that should be performed
        subprocess.call(cmd_to_be_run,shell=True)
        # Delete unnecessary files
        self.__delete_depfiles()
        for f in files_to_be_removed:
            os.remove(f)

    def __get_red_text(self,msg):
        return "\033[31;1m" + msg + "\033[m" # red

    def __get_green_text(self,msg):
        return "\033[32;1m" + msg + "\033[m" # green

#================================
#   Main function
#================================
if __name__ == '__main__':
    tst = Automatic_Tester()
    # Create test configurations
    tst.create_tests()
    # Read the original Makefiles
    tst.read_mkfiles()
    # Read some required files
    tst.read_depfiles()
    # Perform test
    ret = tst.test()
    sys.exit(ret)

