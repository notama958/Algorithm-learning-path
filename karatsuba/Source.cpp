#include "karatsuba.h"

int main() {
	string x = "3141592653589793238462643383279502884197169399375105820974944592";
	string y = "2718281828459045235360287471352662497757247093699959574966967627";
	
	if(x>y)
		cout << karatsuba(x,y) << endl;
	else
		cout << karatsuba(y,x) << endl;
		
	
	return 0;
}