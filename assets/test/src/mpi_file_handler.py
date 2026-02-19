#!/usr/bin/env python3
import socket
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
hostname = socket.gethostname()

# ---------------------------------------------------------
# STEP 1: Rank 0 reads the file and prepares the chunks
# ---------------------------------------------------------
if rank == 0:
    print(f"Rank 0 is reading the file on {hostname}...")
    with open("src/tasks.txt", "r") as f:
        all_lines = f.readlines()
    
    # Strip the hidden newline characters (\n)
    all_lines = [line.strip() for line in all_lines]
    
    # Divide the 8 lines into 4 equal chunks (1 chunk per worker)
    # Result: [['Alpha', 'Bravo'], ['Charlie', 'Delta'], ...]
    chunk_size = len(all_lines) // size
    data_chunks = [all_lines[i : i + chunk_size] for i in range(0, len(all_lines), chunk_size)]
else:
    # Workers must create an empty placeholder variable 
    # so they have somewhere to receive the data!
    data_chunks = None

# ---------------------------------------------------------
# STEP 2: Scatter the chunks over the network
# ---------------------------------------------------------
# Rank 0 deals out the data_chunks. Every rank catches its specific piece in 'my_tasks'.
my_tasks = comm.scatter(data_chunks, root=0)

# ---------------------------------------------------------
# STEP 3: Every Rank processes its specific tasks
# ---------------------------------------------------------
my_results = []
for task in my_tasks:
    # Simulating heavy computation...
    processed_string = f"[{hostname} - Rank {rank}] successfully processed: {task}"
    my_results.append(processed_string)

# ---------------------------------------------------------
# STEP 4: Gather all the processed results back to Rank 0
# ---------------------------------------------------------
all_results = comm.gather(my_results, root=0)

# ---------------------------------------------------------
# STEP 5: Rank 0 writes the final output to a new file
# ---------------------------------------------------------
if rank == 0:
    print("Rank 0 is writing the results to output.txt...")
    with open("src/logs/output.txt", "w") as f:
        for worker_result in all_results:
            for line in worker_result:
                f.write(line + "\n")
    print("Done!")
