FC = gfortran
main=../main

.PHONY: test diff

test1: 
	time $(main) vdw.par test1.cs >& test1_main.out & 
test2:
	time $(main) vdw.par test2.cs >& test2_main.out & 
test3:
	time $(main) vdw.par test3.cs >& test3_main.out & 	
test4:
	time $(main) vdw.par test4.cs >& test4_main.out & 	
test5:
	time $(main) vdw.par test5.cs >& test5_main.out & 	
test6: 
	time $(main) vdw.par 1LYZ_edit_HIS_to_HIP_Q5.cs >& test6_main.out &

diff1: 	
	diff --side-by-side --suppress-common-lines test1.log test1_main.out 
diff2:	
	diff --side-by-side --suppress-common-lines test2.log test2_main.out 
diff3:
	diff --side-by-side --suppress-common-lines test3.log test3_main.out 	
diff4:
	diff --side-by-side --suppress-common-lines test4.log test4_main.out 	
diff5:
	diff --side-by-side --suppress-common-lines test5.log test5_main.out 	
