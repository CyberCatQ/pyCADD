echo
echo 'GPU Info:'
nvidia-smi
echo
read -p 'Enter GPU Number to Use:' GPU
export CUDA_VISIBLE_DEVICES=$GPU
nohup ./pyPMEMD.py > pyPMEMD.out 2>&1 &
echo 'Process Running... You can now safely close the connection.'