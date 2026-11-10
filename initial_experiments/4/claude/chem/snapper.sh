#!/bin/bash
# Copy pr3_r2.npz when the prune log reports target block sizes.
for t in 479 447 415 383 351; do
  until grep -q "\"B\": $t," pr3_r2.log; do sleep 20; done
  sleep 5; cp pr3_r2.npz snap_B$t.npz; echo "snap $t"
done
