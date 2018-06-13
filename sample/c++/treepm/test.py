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
        self.__NPROCS = 4

        # This private variable stores the values of CC, CFLAGS,
        # use_phantom_grape_x86, use_gpu_cuda in the original Makefile.
        self.__mkvars = []

        # This private variable stores the original Makefile.
        self.__mkfile_org = []

    def create_tests(self):
        # Set test configurations
        candidates = []
        # (1) w/ MPI
        CC = ["time mpicc"]
        CXX = ["time mpicxx"]
        CFLAGS = ["-O3 -Wall -DPARTICLE_SIMULATOR_MPI_PARALLEL", \
                  "-O3 -Wall -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp"]
        use_phantom_grape_x86 = ["no", "yes"]
        RUNOPT = ["export OMP_NUM_THREADS={0}; mpirun --mca btl ^openib -np {1} ".format(self.__NTHRDS,self.__NPROCS)]
        listtmp = self.__get_dictlist(CC,CXX,CFLAGS,use_phantom_grape_x86,RUNOPT)
        candidates.extend(listtmp)
        # Remove exceptions
        for i in range(0,len(candidates)):
            item = candidates[i]
            self.__mkvars.append(item) 
        # Check [for debug]
        #for item in self.__mkvars:
        #    print(item)

    def read_mkfile(self):
        # Set filename
        fname = "Makefile"
        # Open file
        fp = open(fname,'r')
        # Read the file
        for line in fp:
            # Delete LF
            #line = line.rstrip()
            # Save the line
            self.__mkfile_org.append(line)
        # Close the file
        fp.close()

    def test(self):
        # Define various constants and declare lists
        pat_set_vars_1st = re.compile(r"^\s*#\s*fdps-autotest-set-vars-1")
        pat_set_vars_2nd = re.compile(r"^\s*#\s*fdps-autotest-set-vars-2")
        pat_run_inst     = re.compile(r"^\s*#\s*fdps-autotest-run")
        param_file_for_test = "param_file_for_test.txt"
        mkfile_for_test  = "Makefile.test"
        logfile_for_test = "stdout-test.txt"
        results = []
        # Perform test
        for mkvars in self.__mkvars:
            # Delete Makefile.test if it exists
            if os.path.isfile(mkfile_for_test):
                os.remove(mkfile_for_test)
            # Create a new Makefile
            fp = open(mkfile_for_test,'w')
            # Write to the file
            for line in self.__mkfile_org:
                # Copy
                fp.write(line)
                # Overwrite the makefile variables if this line
                # corresponds to the insertion point.
                # [1] CC,CXX,CFLAGS
                if (pat_set_vars_1st.match(line)):
                    text = """
                    #---------------------------------------------------
                    # !!! Overwrite the makefile variable for test !!!
                    CC = {CC}
                    CXX = {CXX}
                    CFLAGS = {CFLAGS}
                    #---------------------------------------------------
                    """.format(CC=mkvars["CC"], \
                               CXX=mkvars["CXX"], \
                               CFLAGS=mkvars["CFLAGS"])
                    text = text[1:].rstrip()
                    text = textwrap.dedent(text)
                    text = text + "\n"
                    fp.write(text)
                # [2] use_phantom_grape_x86
                if (pat_set_vars_2nd.match(line)):
                    text = """
                    #---------------------------------------------------
                    # !!! Overwrite the makefile variable for test !!!
                    use_phantom_grape_x86 = {PGFLAG}
                    #---------------------------------------------------
                    """.format(PGFLAG=mkvars["use_phantom_grape_x86"])
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
                    \t{RUNOPT} ./treepm {PARAM_FILE} > {LOGFILE} 2>&1
                    \techo $$? >> {LOGFILE} 2>&1
                    #---------------------------------------------------
                    """.format(RUNOPT=mkvars["RUNOPT"], \
                               PARAM_FILE=param_file_for_test, \
                               LOGFILE=logfile_for_test)
                    text = text[1:].rstrip()
                    text = textwrap.dedent(text)
                    text = text + "\n"
                    fp.write(text)
            # Close the file
            fp.close()
            # Make
            cmd = "make -B -f {0}; echo $?".format(mkfile_for_test)
            outs,errs = subprocess.Popen(cmd, \
                                         stdout=subprocess.PIPE, \
                                         stderr=subprocess.PIPE, \
                                         shell=True).communicate()
            if (outs.rstrip() == "2"): # see https://linuxjm.osdn.jp/html/GNU_make/man1/make.1.html
                msg = self.__get_red_text("failed to make. Stop this test.")
                print("{0}".format(msg))
                print("{0}".format(errs))
                self.__end_processing(results,[mkfile_for_test])
                continue # Move to next test
            # Run
            cmd = "make run -f {0}".format(mkfile_for_test)
            subprocess.call(cmd,shell=True)
            # Obtain a string representing the termination status
            cmd = "tail -n 1 {0}".format(logfile_for_test)
            stat = subprocess.Popen(cmd, \
                                    stdout=subprocess.PIPE, \
                                    shell=True).communicate()[0]
            # Delete all kinds of blank characters from stat
            stat = stat.rstrip()    
            stat = stat.strip()     
            stat = stat.strip("\t")
            # Convert a string to a INT value
            try:
                stat = int(stat)
            except ValueError:
                msg = self.__get_red_text("stat is not a INT value.")
                print("{0}".format(msg))
                self.__end_processing(results,[mkfile_for_test,logfile_for_test])
                continue # Move to next test
            # Store the result
            results.append(stat)
            # Clean
            cmd = "make distclean -f {0}".format(mkfile_for_test)
            subprocess.call(cmd,shell=True)
            # Delete unnecessary files
            os.remove(mkfile_for_test)
            os.remove(logfile_for_test)
            cmd = "rm -rf autotest*"
            subprocess.call(cmd,shell=True)
        # Just in case, we again perform make distclean using the original Makefile
        cmd = "make distclean"
        subprocess.call(cmd,shell=True)
        # Output the summary of the tests
        text = """
        ###################################################
           Summary of tests for sample/c++/treepm
        ###################################################
        """
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
            {TESTNUM}) CC={CC}, CXX={CXX}, CFLAGS={CFLAGS}, use_phantom_grape_x86={PGFLAG}: {STATUS}
            """.format(TESTNUM=testnum, \
                       CC=mkvars["CC"], \
                       CXX=mkvars["CXX"], \
                       CFLAGS=mkvars["CFLAGS"], \
                       PGFLAG=mkvars["use_phantom_grape_x86"], \
                       STATUS=msg)
            text = text[1:].rstrip()
            text = textwrap.dedent(text)
            print("{0}".format(text))
        # Check if all the tests are passed or not.
        return sum(results)

    def __get_dictlist(self,CC,CXX,CFLAGS,use_phantom_grape_x86,RUNOPT):
        listtmp = list(itertools.product(CC, \
                                         CXX, \
                                         CFLAGS, \
                                         use_phantom_grape_x86, \
                                         RUNOPT))
        dictlist = []
        for cc,cxx,cflags,pgflag,runopt in listtmp:
            item = collections.OrderedDict()
            item["CC"] = cc
            item["CXX"] = cxx
            item["CFLAGS"] = cflags
            item["use_phantom_grape_x86"] = pgflag
            item["RUNOPT"] = runopt
            dictlist.append(item)
        return dictlist

    def __end_processing(self,results,files_to_be_removed):
        results.append(1)
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
    # Read the original Makefile
    tst.read_mkfile()
    # Perform test
    ret = tst.test()
    sys.exit(ret)

