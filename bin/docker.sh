#!/bin/bash


set -x


# directory containing whole build process
WORKDIR=/home/jb/workspace

# name of the development image
WORK_IMAGE=flow123d/f123d_docker

get_dev_dir() 
{
    curr_dir=`pwd`
    project_dir="${curr_dir#${WORKDIR}}"
    project_dir="${project_dir#/}"
    project_dir="${project_dir%%/*}"
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
        running_cont=`docker run -itd -v "${WORKDIR}":"${WORKDIR}" flow123d/flow-libs-dev-dbg`        
        
        # setup user and group
        docker exec ${running_cont} addgroup --gid $G_ID $GNAME
        docker exec ${running_cont} adduser --home "$D_HOME" --shell /bin/bash --uid $U_ID --gid $G_ID --disabled-password --system $UNAME
        
        # add git user
        docker exec ${running_cont} git config --global user.email "jbrezmorf@gmail.com"
        docker exec ${running_cont} git config --global user.name "Jan Brezina"
        
        # add git-completion
        curl https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash -o .git-completion.bash
        docker cp .git-completion.bash ${running_cont}:$D_HOME/
        
        # add bashrc, the prompt in particular
        docker cp $WORKDIR/_bashrc_docker ${running_cont}:$D_HOME/.bashrc
                
        # add pmake script
        docker exec ${running_cont} mkdir "$D_HOME/bin"
        docker cp $HOME/bin/pmake ${running_cont}:$D_HOME/bin
        docker exec ${running_cont} chown jb:jb $HOME/bin/pmake
        
        # add ssh keys
        docker exec ${running_cont} mkdir "$D_HOME/.ssh"
        #docker exec ${running_cont} chown jb:jb $HOME/.ssh
        docker cp $HOME/.ssh/id_rsa     ${running_cont}:$D_HOME/.ssh
        docker cp $HOME/.ssh/id_rsa.pub ${running_cont}:$D_HOME/.ssh
        
        
        docker stop ${running_cont}
        docker commit ${running_cont}  $WORK_IMAGE        
    fi    
}


U_ID=`id -u`
G_ID=`id -g`

get_dev_dir

if [ "$1" == "clean" ]
then
    docker stop `docker ps -aq`
    docker rm `docker ps -aq`
    docker rmi $WORK_IMAGE
    exit
elif [ "$1" == "make" ]
then
    shift
    make_work_image
    docker run  --rm -v "${WORKDIR}":"${WORKDIR}" -w "${WORKDIR}/${project_dir}" -u $U_ID:$G_ID $WORK_IMAGE bash -c "$D_HOME/bin/pmake" $@ 
else
    # interactive
    make_work_image
    docker run --rm -it -v "${WORKDIR}":"${WORKDIR}"  -w "${WORKDIR}/${project_dir}" -u $U_ID:$G_ID $WORK_IMAGE bash
fi


