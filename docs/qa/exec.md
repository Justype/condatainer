# Exec Common Errors


## 

            #     Utils.print_note("If .img cannot be mounted, please ensure there is no other process using it.")
            #     Utils.print_note("On a remote server, you can:")
            #     Utils.print_note(Utils.blue_text(f"  fuser -k {last_img}") + " to kill any processes using the img.")
            #     Utils.print_note("On a SLURM system, you can:")
            #     Utils.print_note(Utils.blue_text(f"  squeue -u $USER") + f" to list your jobs. Then use {Utils.blue_text('scancel <JOBID>')} to cancel the job using the img.")
            #     Utils.print_note("After that, run " + Utils.blue_text(f"  e2fsck {last_img}") + " to check the img integrity and repair any issues.")