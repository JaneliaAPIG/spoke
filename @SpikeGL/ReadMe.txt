
Contents.m documentation file is accurate.


Old retired functions
---------------------
FastSettle
GetCurrentSaveFile


Old -> New function replacements
--------------------------------
ConsoleUnhide -> ConsoleShow
IsAcquiring -> IsRunning
DoQueryMatrix -> GetDAQData
GetDir -> EnumRunDir
GetSaveDir -> GetRunDir
GetSaveFile -> GetRunName
SetSaveDir -> SetRunDir
SetSaveFile -> SetRunName
SetSaving -> SetTrgEnable
StartACQ -> StartRun
StopACQ -> StopRun


New functions
-------------
GetAcqChanCounts
SetAOParams
SetAOEnable
SetDigOut


Changed syntax
--------------
GetDAQData now returns two params [mat,headCt]; where headCt is the
zero-based index of the first timepoint in the matrix. This allows
consecutive fetches.


All other functions
-------------------
Same syntax


IMPORTANT
---------
Use new mex files.


