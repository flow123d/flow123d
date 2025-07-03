# Script that mimicks the release build on the Jenkins CI. 
# Can hlep to debug some build errors happening during release build.
set -x

git_branch=`git rev-parse --abbrev-ref HEAD`
build_dir=build-${git_branch}
environment=gnu
release_tag=4.0.0_xyz

if [ ! -d "../flow123d" ]; then
    echo "Not in flow123d source root."
    exit 1
fi

    
if [ -d ${build_dir} ]; then   
    ls .
    echo 
    echo "Recursive remove of '${build_dir}' and 'lib'?"
    sudo rm -rf ${build_dir}
    sudo rm -rf lib
    ls .
fi

config/build/auto_build.sh rel ${environment} ci ${release_tag}
config/build/tar_build_dir.sh ${build_dir}              
config/build/make_packages.sh ${environment} ci ${release_tag}
config/build/make_integration_test_image.sh ${environment} ${release_tag} 
