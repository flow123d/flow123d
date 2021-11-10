#!/bin/bash


set -x

SCRIPT_DIR=${0%/*}

# directory containing whole build process
WORKDIR=/home/jb

# name of the development image
WORK_IMAGE=flow123d/f123d_docker

get_dev_dir() 
{
    curr_dir=`pwd`
    project_dir="${curr_dir#${WORKDIR}}"    # relative to 'workspace'
    project_dir="${project_dir#/}"          
    #project_dir="${project_dir%%/*}"
}

cp_to_docker () {
    source=$1
    target=$D_HOME/$2
    target_file=$D_HOME/$2/${source##*/}
    docker cp $source ${running_cont}:$target
    if ! docker exec ${running_cont} chown $U_ID:$G_ID $target_file
    then    
        docker exec ${running_cont} chown $U_ID:$G_ID $target
    fi        
}

remove_custom_image()
{
    docker stop `docker ps -aq`
    docker rm `docker ps -aq`
    docker rmi $WORK_IMAGE   
}

make_work_image () 
{
    U_ID=`id -u`
    G_ID=`id -g`
    UNAME=`id -un`
    GNAME=`id -gn`
    D_HOME="/home/$UNAME"

    if ! docker images | grep "$WORK_IMAGE" > /dev/null
    then            
        # setup the container
        running_cont=`docker run -itd -v "${WORKDIR}":"${WORKDIR}" flow123d/flow-dev-gnu-dbg`
        
        # setup user and group
        docker exec ${running_cont} addgroup --gid $G_ID $GNAME
        docker exec ${running_cont} adduser --home "$D_HOME" --shell /bin/bash --uid $U_ID --gid $G_ID --disabled-password --system $UNAME
        docker exec ${running_cont} chown $U_ID:$G_ID $D_HOME
        
        # add git user
        docker exec ${running_cont} git config --global user.email "jbrezmorf@gmail.com"
        docker exec ${running_cont} git config --global user.name "Jan Brezina"
        
        # add git-completion
        curl https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash -o .git-completion.bash
        cp_to_docker .git-completion.bash .
        
        # add bashrc, the prompt in particular        
        cp_to_docker $WORKDIR/_bashrc_docker .bashrc
                
        # add pmake script
        docker exec -u $U_ID:$G_ID ${running_cont} mkdir "$D_HOME/bin"
        cp_to_docker $HOME/bin/pmake bin
        
        # add ssh keys
        docker exec -u $U_ID:$G_ID ${running_cont} mkdir "$D_HOME/.ssh"
        #docker exec ${running_cont} chown jb:jb $HOME/.ssh
        cp_to_docker $HOME/.ssh/id_rsa  .ssh
        cp_to_docker $HOME/.ssh/id_rsa.pub  .ssh
                
        docker stop ${running_cont}
        docker commit ${running_cont}  $WORK_IMAGE        
    fi    
}


U_ID=`id -u`
G_ID=`id -g`

get_dev_dir



if [ "$1" == "clean" ]
then
    remove_custom_image
    exit
elif [ "$1" == "update" ]
then    
    remove_custom_image
    # download flow-dev-gnu-dbg image and recreates cusomized container
    # We can not use prebuild images until they are update on regular basis.
    #docker import http://flow.nti.tul.cz/developer-images/flow-libs-dev-dbg.tar.gz    
    make -C $SCRIPT_DIR/docker/dockerfiles flow-dev-gnu-dbg
    make_work_image
elif [ "$1" == "make" ]
then
    shift
    make_work_image
    docker run  --rm -v "${WORKDIR}":"${WORKDIR}" -w "${WORKDIR}/${project_dir}" -u $U_ID:$G_ID $WORK_IMAGE bash "$D_HOME/bin/pmake" "$@" 
elif [ "$1" == "flow123d" ]
then
    shift
    docker run  --rm -v "${WORKDIR}":"${WORKDIR}" -w `pwd` -u $U_ID:$G_ID $WORK_IMAGE bash -c "${SCRIPT_DIR}/flow123d $*" 
else
    # interactive
    make_work_image    
    docker run --rm -it -v "${WORKDIR}":"${WORKDIR}"  -w `pwd` -u $U_ID:$G_ID $WORK_IMAGE bash
fi


