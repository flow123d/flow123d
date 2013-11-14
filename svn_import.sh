#!/bin/bash

# Collection of various functions handy for importing SVN repository to git.
#
# General usage: ./svn_import <function> <params>
#

############################################################3
# Import SVN repo using svn2git tool.
# starting with rev 805 - precursor of version 1.6.0
function import {
  OPT="-v --revision 805 --notags --authors authors --exclude 'JB-deal.ii-richards' --exclude 'branches/JB-richards-anal' --exclude 'PE-IntepolationFactory' --exclude 'third_party/armadillo.' --exclude '.*[.]pos' --exclude 'third_party/boost.*'"
  if [ "$1" == "rebase" ]
  then
    svn2git --rebase $OPT
  else  
    svn2git https://dev.nti.tul.cz/repos/flow123d $OPT
  fi
} 


################################################################
# Call gc; determine size of the repository
function repo_size {
  git gc
  git bundle create tmp.bundle --all
  du -sh tmp.bundle
} 

################################################################
# prints 'HEAD' if the first parameter is in HEAD (?? not sure if we go through all branches)
function is_file_in_head {
    if git ls-tree -r HEAD | grep -q "$1"
    then
        echo "HEAD"
    else
        echo "old"
    fi 
}                  

################################################################
# List all files larger then  $SIZE_LIMIT together with
# filename, blob hash and compressed size
# Output the list into 'bigtosmall.tmp'
function list_large_files {
  SIZE_LIMIT=200000
  git rev-list --objects --all | sort -k 2 > allfileshas.txt

  git verify-pack -v .git/objects/pack/pack-*.idx | grep -v chain | awk -v lim=$SIZE_LIMIT '$3 >lim' | sort -k3nr >bigobjects.txt
  for SHA in `cut -f 1 -d\  < bigobjects.txt`; do
    file=$(grep $SHA allfileshas.txt | awk '{print $2}')
    echo $file $(is_file_in_head $file) $(grep $SHA bigobjects.txt) >> bigtosmall.tmp
  done;
  cat bigtosmall.tmp | column -t > bigtosmall.txt
  rm bigtosmall.tmp
}

################################################################
# Removes file $1 from the history. Filter branch only from the first commit
# to speedup the process.
# 
# Special treatment for the initial commit, since it has no parent. 
function remove_from_history {
  set -x
  export FILE_TO_DELETE=$1
  SHA=`git log --all --oneline -- $FILE_TO_DELETE | awk '{print $1}' | tail -n 1`
  if [ -n "$SHA" ]
  then 
    initial_commit=`git rev-list --oneline --max-parents=0 HEAD | awk '{print $1}'`  
    if [ "$SHA" == "$initial_commit" ]
    then
      range=""
    else
      range="HEAD ^$SHA"
    fi  
    # go through from the first commit of the file, in all branches
    git filter-branch -f --prune-empty --index-filter 'git rm -rf --cached --ignore-unmatch $FILE_TO_DELETE' --tag-name-filter cat -- --all $range
  fi  
}


################################################################
# Usage:
# 1) call list_large_files
# 2) possibly remove or add some lines in bigtosmall.txt
# 3) call clean_history - this removes every file/directory in first column of bigtosmall.txt; that is not in HEAD 
#
# This works, but is rather slow; faster solution could pass just once through the history
# and at every commit delete files from bigtosmall.txt
function clean_history {
        
    # take only those, that are not in HEAD
    cat bigtosmall.txt | grep -v "HEAD" \
    | while read file rest
    do
      remove_from_history $file
    done  
}

##################################################################
# This works very well, but still miss some files.
function clean_history_2 {
  export files=`cat bigtosmall.txt | grep -v "HEAD"` 
  git filter-branch -f --prune-empty --index-filter 'git rm -rf --cached --ignore-unmatch $files' --tag-name-filter cat -- --all
}

#function prune_history {
#  while 
#  list_large_files
#}


function repo_export {
  TARGET=$1
  if [ -n "$TARGET" ]
  then
    rm -rf ${TARGET}
    git clone --no-hardlinks ./.git $TARGET
    cp ./svn_import.sh $TARGET
  fi  
}











################################################################
# Pass through the history and delete big files. Fails since we have to delete the file also in 
# subsequent commits.
#
#
# Pass through all commits (in all branches), at every commit (save the HEAD) delete all files larger then $1. 
# This should work, but it doesn't. Namelly the largest file test_units/mesh/pokus7.msh, is not part of tree at HEAD so
# git rm is not called for it, which is the only difference to:
# git filter-branch ... --index-filter 'git rm ... -rf test_units/mesh/pokus7.msh' -- --all
# ... which works.
#
# Conclusion: use list_large_files or list_large_files_full to select files;
# then delete them through whole history (exclude files present at HEAD)
# possibly can be combined with remove_large_files - which also has some impact on the size (150 -> 120MB)
#
function remove_large_files {
  export SIZE_LIMIT=$1
  echo "SIZE LIMIT: $SIZE_LIMIT"
            # only big files in tree
            # ls-tree fields: <mode> SP <type> SP <object> SP <object size> TAB <file>
              
  git filter-branch  --prune-empty --tree-filter '
        #set -x;  
        head_id=`git rev-parse HEAD`;
        
        TREE_BIG_FILES="`git ls-tree -rl $commit | awk -v lim=$SIZE_LIMIT "\\$4 > lim" | grep -v "/gtest-1.6.0" | grep -v "src/" | awk "{print \\$5}"`";
        git rm -rf --cached --ignore-unmatch $TREE_BIG_FILES;
        
        if [ "$GIT_COMMIT" != "$head_id" ];
        then
            git rm -rf --cached --ignore-unmatch third_party/* >/dev/null;
        fi    
        ' --tag-name-filter cat -- --all
}





############################################################
# Inefficient version of git rev-list ...
function get_blob_refs {
  obj_name="$1"
  shift
  git log "$@" --pretty=format:'%T %h %s' \
  | while read tree commit subject ; do
      if git ls-tree -r $tree | grep -q "$obj_name" ; then
        echo $commit "$subject"
      fi
    done
}

################################################################
# Slow and complex function to get informations provided by list_large_files
# together with all commits wher the file was changed.
# 
# More easily this additional information (if needed) is provided by commad:
# git rev-list --all --oneline -- $file_path
#
function list_large_files_full {
    # Shows you the largest objects in your repo's pack file.
    # (use git gc before)
    #
    # @see http://stubbisms.wordpress.com/2009/07/10/git-script-to-show-largest-pack-objects-and-trim-your-waist-line/
    # @author Antony Stubbs
    
    # set the internal field spereator to line break, so that we can iterate easily over the verify-pack output
    
    # list all objects including their size, sort by size, take top 10
    objects=``

    echo "All sizes are in kB's. The pack column is the size of the object, compressed, inside the pack file."
    
    output="size,pack,SHA,location"

    
    min_size=200000
    # verify-pack fields: sha object_type size compressed 
    
    #BIG_FILES="`git verify-pack -v .git/objects/pack/pack-*.idx | grep -v chain | awk -v lim=$min_size '$3 >$lim' | sort -k3nr`"
    BIG_FILES="`cat xx | grep -v chain `" #| awk -v lim=$min_size '$3 >lim' | sort -k3nr`"
    
    rm -f missing_big
    rm -f size_mismatch
    rm -f found_files
    
    # go through commits
    git log --pretty=format:'%T %H' \
    | while read tree commit remains
    do
        # only big files in tree
        # ls-tree fields: <mode> SP <type> SP <object> SP <object size> TAB <file>
        TREE="`git ls-tree -rl $commit | awk -v lim=$min_size '$4 > lim' | sort -k4nr`"

        # list files changed in current commit
        git show --pretty="format:" --name-only $commit | grep -v '^ *$' > FILES

        # go through intersection (changed big files)
        echo "$TREE" | grep -f FILES \
        | while read a b hash size file rest
        do        
              size=$(( $size / 1024 ))

              # get blob size and compressed size from report of packed repository
              line=`echo "$BIG_FILES" | grep "$hash"`
              fields=( $line )
              if [ "${#fields[@]}" -gt "0" ]
              then

                  # detect if the file is in HEAD
                  if git ls-tree -r HEAD | grep -q $hash
                  then
                      HEAD="HEAD"
                  else
                      HEAD="old"
                  fi  
                  
                  pack_size=$(( ${fields[2]} / 1024))
                  compressed_size=$(( ${fields[3]} / 1024 ))
                  
                  echo "$size,$pack_size,$compressed_size,$file,$HEAD,$commit" >>found_files
#                  echo "$line" >> big_found
              else 
                  echo "$commit $hash $size $file" >>missing_big
              fi   
        done
    done    
    cat found_files | column -t -s ", " | sort -nr > potentialy_big
}


#########################################################################  
 
FUNC=$1
shift
${FUNC} "$@"


#git filter-branch -f --prune-empty --index-filter 'git rm -rf --cached --ignore-unmatch doc/reference_manual/tests_graphics/mel_167s-170_velocity.pdf' --tag-name-filter cat -- 2c56273..
