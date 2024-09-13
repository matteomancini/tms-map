#include <stdio.h>
#include <stdlib.h>


int **convert2d(int *vec, size_t l1, size_t l2);


int seek_seq(int **data, int x_start, int x_end, int q, int y_start, int y_end, int r);


int *lz_complexity(int *svec, size_t l1, size_t l2){

	int c = 1;
	int r = 1;
	int q = 1;
	int k = 1;
	int i = 1;
	int found = 0;
	int a;
	setbuf(stdout, NULL);

	int **ss = convert2d(svec, l1, l2);

	size_t matSize = sizeof(int) * l2;
	int *carray = (int*)(malloc(matSize));

	while(r <= l2){
		if(q == r) {
			a = i+k-1;
		} else {
			a = l1;
		}
		
		found = seek_seq(ss, 0, a, q-1, i, i+k, r-1);

		if(found){
			k += 1;
			if (i + k > l1){
				carray[r-1] = c;
				r += 1;
				printf("%d\n", r);
				i = 0;
				q = r - 1;
				k = 1;
			}
		} else {
			q -= 1;
			if(q < 1){
				c += 1;
				i += k;
				if(i + 1 > l1){
					carray[r-1] = c;
					r += 1;
					printf("%d\n", r);
					i = 0;
					q = r - 1;
					k = 1;
				} else {
					q = r;
					k = 1;
				}
			}
		}
	}
	carray[r-2] += 1;
	return(carray);
}


int **convert2d(int *vec, size_t l1, size_t l2){

	int **mat;

	mat = (int **)malloc(sizeof(int*)*l1);

	for (int i = 0; i < l1; i++) {
		mat[i] = &vec[i*l2];
	}

	return mat;
}


int seek_seq(int **data, int x_start, int x_end, int q, int y_start, int y_end, int r){

	int found = 0;
	for (int i = x_start; (i < x_end - y_end + y_start + 1) && !found; i++){
		for (int j = 0; j < y_end - y_start; j++){
			if(data[j+i][q] != data[y_start + j][r]){
				break;
			}
			if(j == y_end - y_start - 1){
				found = 1;
			}
		}
	}
	return(found);
}

