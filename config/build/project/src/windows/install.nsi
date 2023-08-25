# Flow123d nsi build script for Windows


#----------------------------------------------------------
# GENERAL
#----------------------------------------------------------
  # Name and file
  # Read version information from file.
  !searchparse /file "version" '' VERSION
  !searchparse /file "imagename" '' IMAGE ''

  # Name and file
  Name "Flow123d ${VERSION}"
  Caption "Flow123d ${VERSION} Setup"
  OutFile "flow123d_${VERSION}_windows_install.exe"


  # Admin access
  # RequestExecutionLevel admin

  # Change the default Modern UI icons
  # !define MUI_ICON "${NSISDIR}\Contrib\Graphics\Icons\modern-install-blue.ico"
  # !define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\classic-uninstall.ico"

  # Default installation folder
  InstallDir "$LOCALAPPDATA\Flow123d\v${VERSION}"
  !define UninstName "Uninstall"

  # Maximum compression
  # https://nsis.sourceforge.io/Reference/SetCompressor
  # In our case the final exe:
  #   zlib (default): 58 MB in 20 sec
  #   bzip2:          43 MB in 1 min
  #   lzma:           28 MB in 5 min
  SetCompressor /SOLID zlib

#----------------------------------------------------------
# Interface Settings
#----------------------------------------------------------
  # installation only for current user

  # Include the tools we use.
  !include MUI2.nsh
  !include LogicLib.nsh
  !include WinMessages.nsh
  !include x64.nsh
  !include FileFunc.nsh
  !include nsis\EnvVarUpdate.nsh
  !include nsis\dumplog.nsi
  !include nsis\docker_install.nsi
  !insertmacro GetDrives
  # !include MultiUser.nsh

  !define MUI_ICON "win\logo.ico"
  !define MULTIUSER_EXECUTIONLEVEL Standard
  !define MULTIUSER_INSTALLMODE_INSTDIR   "Flow123d"
  !define MUI_FINISHPAGE_NOAUTOCLOSE
  !define MUI_CUSTOMFUNCTION_ABORT CustomOnUserAbort

  !define MUI_WELCOMEFINISHPAGE_BITMAP    "win\back.bmp" # 170x320
  !define MUI_UNWELCOMEFINISHPAGE_BITMAP  "win\back.bmp" # 170x320

#----------------------------------------------------------
# PAGES
#----------------------------------------------------------
  # Define the different pages.
  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE     "LICENSE.txt"
  #!insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  ; custom Docker page
  Page custom DockerPage DockerPageLeave
  ; for enabling cancel button during installation
  !define MUI_PAGE_CUSTOMFUNCTION_SHOW InstFilesShow
  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH
  # uninstall
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  !insertmacro MUI_LANGUAGE         "English"

#----------------------------------------------------------
# Functions
#----------------------------------------------------------
Function InstFilesShow
  ; see https://nsis.sourceforge.io/InstFiles_Cancel_-_Allowing_a_user_to_cancel_installation_during_InstFiles
  GetDlgItem $0 $HWNDPARENT 2
  EnableWindow $0 1
FunctionEnd

Function COPY_FILES
  SetOutPath $INSTDIR
  File "version"
  File "license.txt"
  File "readme.txt"

  File /r "tests"

  SetOutPath "$INSTDIR\doc"
  File /r "htmldoc"
  File "flow123d_${VERSION}_doc.pdf"

  # install fterm.bat
  CreateDirectory "$INSTDIR\bin"
  SetOutPath "$INSTDIR\bin"

  # fterm.bat serves as an interactive terminal for (mostly) debugging purposes
  # or for more experience users
  FileOpen $0 "fterm.bat" w
    FileWrite $0 '@echo off$\r$\n'
    FileWrite $0 'SET cdir=\%CD:~0,1%\%CD:~3,256%$\r$\n'
    FileWrite $0 'SET L=%CD:~0,1%$\r$\n'
    FileWrite $0 'docker run -ti --rm -v "%L%:\:/%L%/" -v "c:\:/c/" -w "%cdir:\=/%" ${IMAGE} %*$\r$\n'
    FileWrite $0 'pause$\r$\n'
  FileClose $0

  # flow123d.bat is wrapper which directly calls a binary inside a docker image
  # flow123d.bat and flow123d (inside docker) are linked together, meaning when flow123d.bat is over
  # it means the docker flow123d inside fihisned and container exited.
  FileOpen $0 "flow123d.bat" w
    FileWrite $0 '@echo off$\r$\n'
    FileWrite $0 'SET cdir=\%CD:~0,1%\%CD:~3,256%$\r$\n'
    FileWrite $0 'SET L=%CD:~0,1%$\r$\n'
    FileWrite $0 'docker run -ti --rm -v "%L%:\:/%L%/" -v "c:\:/c/" -w "%cdir:\=/%" ${IMAGE} flow123d %*$\r$\n'
  FileClose $0

  # In some cases when executing flow123d from other process we need to start docker without terminal in non-interactive mode.
  # flow-noterm.bat serves this purpose:
  FileOpen $0 "flow-noterm.bat" w
    FileWrite $0 '@echo off$\r$\n'
    FileWrite $0 'SET cdir=\%CD:~0,1%\%CD:~3,256%$\r$\n'
    FileWrite $0 'SET L=%CD:~0,1%$\r$\n'
    FileWrite $0 'docker run --rm -v "%L%:\:/%L%/" -v "c:\:/c/" -w "%cdir:\=/%" ${IMAGE} flow123d %*$\r$\n'
  FileClose $0

  # fterm and flow123d with version
  CopyFiles "fterm.bat"     "fterm-${VERSION}.bat"
  CopyFiles "flow123d.bat"  "flow123d-${VERSION}.bat"
  CopyFiles "flow123d.bat"  "flow.bat"
  CopyFiles "flow123d.bat"  "flow-${VERSION}.bat"  
  CopyFiles "flow-noterm.bat"  "flow-noterm-${VERSION}.bat"

FunctionEnd

Function START_DOCKER_DEAMON
  # try to run docker ps to determine whther deamon is running
  DetailPrint "Docker installed in $DOCKER_PATH"
  nsExec::Exec 'docker ps'
  Pop $0

  # if is running continue
  ${If} $0 == "0"
    DetailPrint "Docker deamon running"
  # otherwise lanch the deamon
  ${Else}
    DetailPrint "Starting Docker deamon, this'll take about a minute for the first time."
    Exec '$DOCKER_EXE'
    Sleep 10000
    MessageBox MB_OK "Docker Desktop is now starting, wait until it's ready before pressing OK$\r$\nYou should see a notification at the bottom right."
  ${EndIf}
FunctionEnd


Function PULL_IMAGE
  # pull image
  DetailPrint "Pulling Docker image: ${IMAGE}"
  MessageBox MB_OK "Docker will now pull the Flow123d image from internet.$\r$\nThis might take some time depending on your connection."
  nsExec::Exec 'docker pull ${IMAGE}'
FunctionEnd


Function MAKE_DOCKERD
    # install dockerd.bat
    SetOutPath "$INSTDIR\bin"
    FileOpen $0 "dockerd.bat" w
      FileWrite $0 '@echo off$\r$\n'
      FileWrite $0 'Start "" "$DOCKER_EXE$\r$\n"'
    FileClose $0
FunctionEnd

Function REGISTER_APP

  # Update path
  ${EnvVarUpdate}  $0 "PATH" "P" "HKCU" "$INSTDIR\bin"
  # Update global FLOW123D env variable version list
  ${EnvVarUpdate}  $0 "FLOW123D" "A" "HKCU" "${VERSION}"

  #----------------------------------------------------------

  # Write the installation path into the registry
  WriteRegStr HKCU "SOFTWARE\Flow123d-${VERSION}" "Install_Dir" "$INSTDIR"
  # Write the uninstall keys for Windows
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\Flow123d-${VERSION}" "DisplayName" "Flow123d-${VERSION}"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\Flow123d-${VERSION}" "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\Flow123d-${VERSION}" "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\Flow123d-${VERSION}" "NoRepair" 1
FunctionEnd

Function DumpLogToFile
  StrCpy $0 "$INSTDIR\install.log"
	Push $0
  Call DumpLog
FunctionEnd

Function CustomOnUserAbort
   MessageBox MB_YESNO "Abort install?" IDYES NoCancelAbort
     Call DumpLogToFile
     Abort
   NoCancelAbort:
FunctionEnd

Function .onInstFailed
  MessageBox MB_OK "Installation failed."
  Call DumpLogToFile
FunctionEnd

#----------------------------------------------------------
# Section
#----------------------------------------------------------
Section

  SetDetailsView show

  RMDir /r $INSTDIR
  Call REMOVE_OLD
  ; Call CHECK_DOCKER
  Call COPY_FILES
  Call START_DOCKER_DEAMON
  Call PULL_IMAGE
  Call MAKE_DOCKERD
  Call REGISTER_APP

  CreateShortcut "$DESKTOP\Flow-${VERSION}.lnk" "$INSTDIR\bin\fterm.bat" "" "$INSTDIR\win\logo-64.ico"

  #----------------------------------------------------------

  WriteUninstaller "uninstall.exe"

  Call DumpLogToFile

SectionEnd



# Uninstaller
Section "Uninstall"

  SetDetailsView show

  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\Flow123d-${VERSION}"
  DeleteRegKey HKCU "SOFTWARE\Flow123d-${VERSION}"
  Delete "$DESKTOP\Flow-${VERSION}.lnk"

  ${un.EnvVarUpdate} $0 "PATH"     "R" "HKCU" "$INSTDIR\bin"
  ${un.EnvVarUpdate} $0 "FLOW123D" "R" "HKCU" "${VERSION}"

  # Remove docker image
  ExecWait 'docker rmi -f ${IMAGE}'
  # Remove Flow123d installation directory and all files within.
  RMDir /r "$INSTDIR"

SectionEnd
