#! /bin/bash	
# Use the bash shell to run the NGH sample
echo
echo "################################################################################"
echo "##########                            FLOW                           ## Start ##"
echo "################################################################################"
echo

# /home/dalibor/NetBeansProjects/Flow123D/1.6.5/bin/flow123d   -s   ./04_flow/flow.ini
 /home/dalibor/bin/flow_1.6.6_linux_64/bin/flow123d   -s   ./04_flow/flow.ini

echo
echo "################################################################################"
echo "##########                            FLOW                           ##  End  ##"
echo "################################################################################"
echo

read -p "Press any key to continue"

