#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=30GB
#SBATCH --gpus nvidia_h100_80gb_hbm3_1g.10gb

#CNT --base-image src/pytorch.sqf

echo CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES
echo NVIDIA_VISIBLE_DEVICES $NVIDIA_VISIBLE_DEVICES
echo NVIDIA_DRIVER_CAPABILITIES $NVIDIA_DRIVER_CAPABILITIES

nvidia-smi
echo "Running PyTorch GPU test..."
python3 -c "import torch; print('PyTorch version:', torch.__version__); print('CUDA available:', torch.cuda.is_available())"
