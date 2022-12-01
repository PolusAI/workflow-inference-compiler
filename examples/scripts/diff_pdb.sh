# cut out the atom number column first, so renumberings (i.e. due to inserted/deleted atoms) don't cause spurious diffs
diff -u <(sed -r 's/(\s+)?\S+//2' $1) <(sed -r 's/(\s+)?\S+//2' $2) || :
# NOTE: Unbelievably, diff returns 0 if no differences were found and 1 if differences were found.
# Use || : to always return 0 so we don't make the CWL runner think the script failed.