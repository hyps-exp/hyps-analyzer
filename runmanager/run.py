#!/opt/python-3.5/bin/python

#____________________________________________________

__author__  = 'Y.Nakada <nakada@km.phys.sci.osaka-u.ac.jp>'
__version__ = '2.0'
__date__    = '15 April 2018'

#____________________________________________________

import os
import sys
import time
import signal

sys.path.append( os.path.dirname( os.path.abspath( sys.argv[0] ) )
                 + '/module' )

import utility
import RunlistManager
import JobManager

#____________________________________________________

SLEEP_TIME = 3

#____________________________________________________

def handler( num, frame ) :

    print( 'KeyboardInterrupt' )

    print( 'Terminating processes...' )
    for index, item in enumerate( joblist ) :
        item[0].killAll()
        joblist[index][1] = 'terminated'

    utility.updateJobStat( joblist, fjobstat )

    time.sleep( 1 )             # waiting for bsub log files are generated

    print( 'Deleting files...' )
    for item in joblist :
        # item[0].clear()
        item[0].clearAll()

    sys.exit( 0 )

#____________________________________________________

def decodeTime( second ) :

    second = int( second )

    hour    = second // 3600
    second -= hour * 3600
    minute  = second // 60
    second -= minute * 60

    return hour, minute, second

#____________________________________________________

argvs = sys.argv
argc = len( argvs )

if argc != 2 :
    print( 'USAGE: %s [ file ]' % argvs[0] )
    sys.exit( 0 )

if not os.path.exists( argvs[1] ) :
    utility.ExitFailure( 'No such file: %s' % argvs[1] )

frunlist = argvs[1]

#____________________________________________________

djobstat = os.path.dirname( os.path.abspath( sys.argv[0] ) ) + '/stat'
if not os.path.exists( djobstat ) :
    os.mkdir( djobstat )

fjobstat = djobstat + '/'\
           + os.path.splitext( os.path.basename( frunlist ) )[0]\
           + '.stat'

#____________________________________________________

runMan  = RunlistManager.RunlistManager( frunlist )
runtag  = runMan.getTag()
runlist = runMan.getRunlist()

joblist = list()

for run in runlist :
    jobMan = JobManager.JobManager( runtag, *run )
    joblist.append( [ jobMan, 'init', decodeTime( jobMan.getDiffTime() ) ] )

signal.signal( signal.SIGINT, handler )

utility.updateJobStat( joblist, fjobstat )

njobs = len( joblist )
fl_done = 0

while not fl_done  == njobs :
    for index, ( job, stat, ptime ) in enumerate( joblist ) :
        if stat == 'init' :
            job.execute()
            joblist[index][1] = 'running(0/%d)' % job.getNSegs()
            joblist[index][2] = decodeTime( job.getDiffTime() )
        elif stat[:7] == 'running' :
            job.updateJobStat()
            if job.isExecuted() :
                job.mergeFOut()
                joblist[index][1] = 'merging'
            else :
                numer = job.getProgress()
                denom = job.getNSegs()
                joblist[index][1] = 'running(%d/%d)' % ( numer, denom )
            joblist[index][2] = decodeTime( job.getDiffTime() )
        elif stat == 'merging' :
            fstat = job.getFinalStatus()
            if fstat is None :
                continue
            elif fstat is True :
                # job.clear()
                job.clearAll()
                joblist[index][1] = 'done'
            elif fstat is False :
                job.clear()
                joblist[index][1] = 'error'
            else :
                job.clear()
                joblist[index][1] = 'unknown'
            fl_done += 1
            joblist[index][2] = decodeTime( job.getDiffTime() )
        elif stat == 'done' :
            continue
        elif stat == 'error' :
            continue
        else :
            joblist[index][1] = 'unknown'
            joblist[index][2] = decodeTime( job.getDiffTime() )

    utility.updateJobStat( joblist, fjobstat )
    time.sleep( SLEEP_TIME )
