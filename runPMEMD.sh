echo
echo 'GPU Info:'
nvidia-smi
echo
read -p 'Enter GPU Number to Use:' GPU
export CUDA_VISIBLE_DEVICES=$GPU
nohup pyPMEMD.py > pyMEMD.out 2>&1 &