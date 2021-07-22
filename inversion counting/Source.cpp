#include "inversion.h"

int main()
{
	vector<int>vt = fileIn();
	//loop(vt);
	cout << "Inversion counting: " << inversion_counting(vt, 0, vt.size() - 1) << endl;
	return 0;
}