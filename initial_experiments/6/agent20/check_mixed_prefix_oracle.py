"""Use the complete native small-oracle contract with the mixed reference."""
import argparse
import check_prefix_counter_oracle as check
import mixed_prefix_oracle as oracle
import mixed_prefix_reference as reference
check.oracle=oracle;check.reference=reference
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=2048);a=p.parse_args();check.test(pairs=2,through=2,shots=a.shots)
