#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

// file input => vector<int>
vector<int> fileIn(string file = "data.txt")
{
	// input file should not contain empty line
	ifstream fin;
	fin.open(file);
	vector<int>result;
	if (fin.is_open())
	{
		int x;
		while (fin >> x)
		{
			result.push_back(x);
		}
	}

	return result;
}
void loop(vector<int> vt)
{
	for (int i : vt)
		cout << i << ", ";
	cout << endl;
}

// merge two sorted arrays
long long  merge(vector<int>&arr, int left, int mid, int right)
{

	long long  count = 0;
	int subArrLeft = mid - left + 1;
	int subArrRight = right - mid;
	// create 2 sub vector arrs
	vector<int> left_arr;
	for (int i = 0; i < subArrLeft; i++)
	{
		left_arr.push_back(arr[left+i]);
	}
	/*loop(left_arr);
	cout << "----" << endl;*/
	vector<int> right_arr;
	for (int i = mid+1; i <= right; i++)
	{
		right_arr.push_back(arr[i]);
	}
	//loop(left_arr);
	//loop(right_arr);
	int i = 0, j = 0;
	int k = left;
	vector<int> temp(arr.size());
	while (i < subArrLeft && j < subArrRight)
	{
		if (left_arr[i] <= right_arr[j])
		{
			temp[k] = left_arr[i];
			i++;
		}
		else {
			temp[k] = right_arr[j];
			j++;
			// eg input (1,20,6,4,5)
			// first count increment +1
			// (1,20) and 6
			// 1,6	20
			// second count increment +2
			// (1,6,20) (4,5)
			// 1,4
			// third/final count increment +2
			// (1,6,20) (4,5)
			// 1,4,5	6 20
			//cout << "Herer:"<<subArrLeft-i<< endl;
			count += (subArrLeft - i);
		}
		k++;
	}
	// merge the remaining back to arr
	while (i < subArrLeft)
	{
		temp[k++] = left_arr[i++];
	}
	while (j < subArrRight)
	{
		temp[k++] = right_arr[j++];
	}
	// merge sorted els to arr
	for (int i = left; i <= right; i++)
	{
		arr[i] = temp[i];
	}
	/*loop(temp);
	loop(arr);
	cout << count << endl;*/
	return count;
}
// divide and merge and sum up inversion
// right and left are index
long long inversion_counting(vector<int> &arr, int left, int right)
{
	int mid= 0;
	long long count = 0;
	if (right > left)
	{
		mid = (right + left) / 2;
		count += inversion_counting(arr, left, mid);
		count += inversion_counting(arr, mid + 1, right);
		// merge
		count += merge(arr, left, mid, right);
	}
	return count;
}



