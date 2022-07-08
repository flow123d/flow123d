# see https://docs.docker.com/desktop/windows/install/
$url = "https://desktop.docker.com/win/main/amd64/Docker Desktop Installer.exe"
$output = "dfw.exe"
Start-BitsTransfer -Source $url -Destination $output
