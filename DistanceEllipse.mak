# Microsoft Developer Studio Generated NMAKE File, Format Version 4.20

# ** DO NOT EDIT **



# TARGTYPE "Win32 (x86) Console Application" 0x0103



!IF "$(CFG)" == ""

CFG=DistanceEllipse - Win32 Debug

!MESSAGE No configuration specified.  Defaulting to DistanceEllipse - Win32\

 Debug.

!ENDIF 



!IF "$(CFG)" != "DistanceEllipse - Win32 Release" && "$(CFG)" !=\

 "DistanceEllipse - Win32 Debug"

!MESSAGE Invalid configuration "$(CFG)" specified.

!MESSAGE You can specify a configuration when running NMAKE on this makefile

!MESSAGE by defining the macro CFG on the command line.  For example:

!MESSAGE 

!MESSAGE NMAKE /f "DistanceEllipse.mak" CFG="DistanceEllipse - Win32 Debug"

!MESSAGE 

!MESSAGE Possible choices for configuration are:

!MESSAGE 

!MESSAGE "DistanceEllipse - Win32 Release" (based on\

 "Win32 (x86) Console Application")

!MESSAGE "DistanceEllipse - Win32 Debug" (based on\

 "Win32 (x86) Console Application")

!MESSAGE 

!ERROR An invalid configuration is specified.

!ENDIF 



!IF "$(OS)" == "Windows_NT"

NULL=

!ELSE 

NULL=nul

!ENDIF 

################################################################################

# Begin Project

# PROP Target_Last_Scanned "DistanceEllipse - Win32 Debug"

RSC=rc.exe

CPP=cl.exe



!IF  "$(CFG)" == "DistanceEllipse - Win32 Release"



# PROP BASE Use_MFC 0

# PROP BASE Use_Debug_Libraries 0

# PROP BASE Output_Dir "Release"

# PROP BASE Intermediate_Dir "Release"

# PROP BASE Target_Dir ""

# PROP Use_MFC 0

# PROP Use_Debug_Libraries 0

# PROP Output_Dir "Release"

# PROP Intermediate_Dir "Release"

# PROP Target_Dir ""

OUTDIR=.\Release

INTDIR=.\Release



ALL : "$(OUTDIR)\DistanceEllipse.exe"



CLEAN : 

	-@erase "$(INTDIR)\DrawEllipse.obj"

	-@erase "$(INTDIR)\EstimateNear.obj"

	-@erase "$(INTDIR)\FindNear.obj"

	-@erase "$(OUTDIR)\DistanceEllipse.exe"



"$(OUTDIR)" :

    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"



# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /YX /c

# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /YX /c

CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE"\

 /Fp"$(INTDIR)/DistanceEllipse.pch" /YX /Fo"$(INTDIR)/" /c 

CPP_OBJS=.\Release/

CPP_SBRS=.\.

# ADD BASE RSC /l 0x409 /d "NDEBUG"

# ADD RSC /l 0x409 /d "NDEBUG"

BSC32=bscmake.exe

# ADD BASE BSC32 /nologo

# ADD BSC32 /nologo

BSC32_FLAGS=/nologo /o"$(OUTDIR)/DistanceEllipse.bsc" 

BSC32_SBRS= \

	

LINK32=link.exe

# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\

 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\

 odbccp32.lib /nologo /subsystem:console /incremental:no\

 /pdb:"$(OUTDIR)/DistanceEllipse.pdb" /machine:I386\

 /out:"$(OUTDIR)/DistanceEllipse.exe" 

LINK32_OBJS= \

	"$(INTDIR)\DrawEllipse.obj" \

	"$(INTDIR)\EstimateNear.obj" \

	"$(INTDIR)\FindNear.obj"



"$(OUTDIR)\DistanceEllipse.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)

    $(LINK32) @<<

  $(LINK32_FLAGS) $(LINK32_OBJS)

<<



!ELSEIF  "$(CFG)" == "DistanceEllipse - Win32 Debug"



# PROP BASE Use_MFC 0

# PROP BASE Use_Debug_Libraries 1

# PROP BASE Output_Dir "Debug"

# PROP BASE Intermediate_Dir "Debug"

# PROP BASE Target_Dir ""

# PROP Use_MFC 0

# PROP Use_Debug_Libraries 1

# PROP Output_Dir "Debug"

# PROP Intermediate_Dir "Debug"

# PROP Target_Dir ""

OUTDIR=.\Debug

INTDIR=.\Debug



ALL : "$(OUTDIR)\DistanceEllipse.exe"



CLEAN : 

	-@erase "$(INTDIR)\DrawEllipse.obj"

	-@erase "$(INTDIR)\EstimateNear.obj"

	-@erase "$(INTDIR)\FindNear.obj"

	-@erase "$(INTDIR)\vc40.idb"

	-@erase "$(INTDIR)\vc40.pdb"

	-@erase "$(OUTDIR)\DistanceEllipse.exe"

	-@erase "$(OUTDIR)\DistanceEllipse.ilk"

	-@erase "$(OUTDIR)\DistanceEllipse.pdb"



"$(OUTDIR)" :

    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"



# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /YX /c

# ADD CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /YX /c

CPP_PROJ=/nologo /MLd /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE"\

 /Fp"$(INTDIR)/DistanceEllipse.pch" /YX /Fo"$(INTDIR)/" /Fd"$(INTDIR)/" /c 

CPP_OBJS=.\Debug/

CPP_SBRS=.\.

# ADD BASE RSC /l 0x409 /d "_DEBUG"

# ADD RSC /l 0x409 /d "_DEBUG"

BSC32=bscmake.exe

# ADD BASE BSC32 /nologo

# ADD BSC32 /nologo

BSC32_FLAGS=/nologo /o"$(OUTDIR)/DistanceEllipse.bsc" 

BSC32_SBRS= \

	

LINK32=link.exe

# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386

# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386

LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib\

 advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib\

 odbccp32.lib /nologo /subsystem:console /incremental:yes\

 /pdb:"$(OUTDIR)/DistanceEllipse.pdb" /debug /machine:I386\

 /out:"$(OUTDIR)/DistanceEllipse.exe" 

LINK32_OBJS= \

	"$(INTDIR)\DrawEllipse.obj" \

	"$(INTDIR)\EstimateNear.obj" \

	"$(INTDIR)\FindNear.obj"



"$(OUTDIR)\DistanceEllipse.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)

    $(LINK32) @<<

  $(LINK32_FLAGS) $(LINK32_OBJS)

<<



!ENDIF 



.c{$(CPP_OBJS)}.obj:

   $(CPP) $(CPP_PROJ) $<  



.cpp{$(CPP_OBJS)}.obj:

   $(CPP) $(CPP_PROJ) $<  



.cxx{$(CPP_OBJS)}.obj:

   $(CPP) $(CPP_PROJ) $<  



.c{$(CPP_SBRS)}.sbr:

   $(CPP) $(CPP_PROJ) $<  



.cpp{$(CPP_SBRS)}.sbr:

   $(CPP) $(CPP_PROJ) $<  



.cxx{$(CPP_SBRS)}.sbr:

   $(CPP) $(CPP_PROJ) $<  



################################################################################

# Begin Target



# Name "DistanceEllipse - Win32 Release"

# Name "DistanceEllipse - Win32 Debug"



!IF  "$(CFG)" == "DistanceEllipse - Win32 Release"



!ELSEIF  "$(CFG)" == "DistanceEllipse - Win32 Debug"



!ENDIF 



################################################################################

# Begin Source File



SOURCE=.\DrawEllipse.cpp

DEP_CPP_DRAWE=\

	".\EstimateNear.h"\

	".\FindNear.h"\

	



"$(INTDIR)\DrawEllipse.obj" : $(SOURCE) $(DEP_CPP_DRAWE) "$(INTDIR)"





# End Source File

################################################################################

# Begin Source File



SOURCE=.\FindNear.cpp



"$(INTDIR)\FindNear.obj" : $(SOURCE) "$(INTDIR)"





# End Source File

################################################################################

# Begin Source File



SOURCE=.\EstimateNear.cpp

DEP_CPP_ESTIM=\

	".\EstimateNear.h"\

	



"$(INTDIR)\EstimateNear.obj" : $(SOURCE) $(DEP_CPP_ESTIM) "$(INTDIR)"





# End Source File

# End Target

# End Project

################################################################################

