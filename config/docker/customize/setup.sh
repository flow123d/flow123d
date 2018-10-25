#!/bin/bash
# author: Jan Hybs
#
# this is template file for following docker images
#  - flow-libs-dev-dbg
#  - flow-libs-dev-rel
#
# purpose of this file is to create more easy to use
# environment in docker images
# variables surrounded with @ will be replaced later, available variables are:
#     uname     - string username
#     gid       - group id
#     uid       - user id
#     git_email - result from `git config --global user.email`
#     uname - result from `git config --global user.name`
#
# this script will be executed inside running docker container
# right now you are a root with unlimited privileges



# CREATE USER AND GROUP
# ------------------------------------------------------------------------------
echo "Running customization script..."
echo "-------------------------------"
echo 'group: @uname@(@gid@)'
echo 'user:  @uname@(@uid@)'

addgroup  --gid @gid@ --force-badname @uname@
adduser   --home /home/@uname@ --shell /bin/bash \
          --uid @uid@ --gid @gid@ \
          --disabled-password --system --force-badname @uname@


# BUILDER COMMANDS
# ------------------------------------------------------------------------------
# create folder where user will have access to
mkdir -p /opt/flow123d/flow123d
chown -R @uid@:@gid@ /opt/flow123d/

# allow sudo for user
cat >> /etc/sudoers  << EOL
@uname@ ALL=NOPASSWD: ALL
EOL

# edit main bash.bashrc file
cat >> /etc/bash.bashrc << EOL
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi
# shortcuts
alias ls='ls --color=auto'
alias grep='grep --color=auto'
alias fgrep='fgrep --color=auto'
alias egrep='egrep --color=auto'

function _flow123d() {
  local cur=\${COMP_WORDS[COMP_CWORD]}
  COMPREPLY=(\$(compgen -df -W "-s -i -o -l \
        --solve --input_dir --output_dir --log --version \
        --no_log --no_signal_handler --no_profiler --help \
        --input_format --petsc_redirect --yaml_balance" -- \$cur))
  return 0
}

function _runtest() {
  local cur=\${COMP_WORDS[COMP_CWORD]}
  COMPREPLY=(\$(compgen -df -W "--keep-going \
        --batch --cpu --limit-time --limit-memory \
        --no-clean --no-compare --death-test --help \
        --list" -- \$cur))
  return 0
}
complete -o nospace -F _flow123d flow123d
complete -o nospace -F _runtest runtest

EOL


# edit .bashrc ($PS1 variable) so the version is visible
image_tag="@image_tag@"
if [[ -z "${image_tag}" ]]; then
    image_tag="@docker_image@"
fi
image_version="@image_version@"
if [[ -z "${image_version}" ]]; then
    image_version="rel"
fi
image_version="$image_version-"
image_id=${image_version}${image_tag}

# Color definition
txtblk='\e[0;30m' # Black - Regular
txtred='\e[0;31m' # Red
txtgrn='\e[0;32m' # Green
txtylw='\e[0;33m' # Yellow
txtblu='\e[0;34m' # Blue
txtpur='\e[0;35m' # Purple
txtcyn='\e[0;36m' # Cyan
txtwht='\e[0;37m' # White
bldblk='\e[1;30m' # Black - Bold
bldred='\e[1;31m' # Red
bldgrn='\e[1;32m' # Green
bldylw='\e[1;33m' # Yellow
bldblu='\e[1;34m' # Blue
bldpur='\e[1;35m' # Purple
bldcyn='\e[1;36m' # Cyan
bldwht='\e[1;37m' # White
txtrst='\e[0m'    # Text Reset


LIGHT_THEME="@light_theme@"
if [[ "${LIGHT_THEME}" -eq 1 ]]; then
  PS1_DOCKER="${bldgrn}\u${txtrst}@${bldgrn}${image_id}${bldylw} \w ${txtrst}"
else
  PS1_DOCKER="${bldpur}\u${txtrst}@${bldpur}${image_id}${bldblu} \w ${txtrst}"
fi

cat >> /etc/bash.bashrc << EOL
export PS1="$PS1_DOCKER"
# clear the terminal
printf '\033[2J'
echo " ___ _            _ ___ ____    _  "
echo "| __| |_____ __ _/ |_  )__ / __| | "
echo "| _|| / _ \ V  V / |/ / |_ \/ _  | "
echo "|_| |_\___/\_/\_/|_/___|___/\__,_| "
echo "                          $image_id"
EOL

# copy git configuration
cp -r /tmp/.gitconfig /home/@uname@/

# copy ssh keys
cp -r /tmp/.ssh /home/@uname@/

# chown files
chown -R @uid@:@gid@ /home/@uname@
