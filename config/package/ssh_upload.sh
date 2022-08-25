source_dir=$1   # local source dir which contents will be uploaded
target_dir=$2   # target image on the remote, will be created
ssh_remote=$3   # user@host
ssh_pass=$4     # password

set -x
relver="${{ steps.vars.outputs.relver }}"
tar -cvf package.tar ${source_dir}

pass_cmd=sshpass -p ${ssh_pass}
opt=-o PreferredAuthentications=password -o PubkeyAuthentication=no -o StrictHostKeyChecking=no

# Make directory with parents (-p also does not complain if the directory exists)
${pass_cmd} ssh ${opt} ${ssh_remote} mkdir -p ${target_dir}

# copy the tar archive
${pass_cmd} scp -vr ${opt} package.tar ${ssh_remote}:${target_dir}/package.tar

# unter archive
${pass_cmd} ssh ${opt} ${ssh_remote} \
tar -xf ${target_dir}/package.tar -C ${target_dir} --strip-components=1

# remove archive
${pass_cmd} ssh ${opt} ${ssh_remote} rm ${target_dir}/package.tar
