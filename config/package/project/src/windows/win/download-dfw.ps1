$url = "http://download.docker.com/win/stable/Docker for Windows Installer.exe"
$output = "dfw.exe"
Start-BitsTransfer -Source $url -Destination $output
