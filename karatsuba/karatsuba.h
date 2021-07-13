#pragma once
#include <iostream>
#include <string>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
string checkTrailingZeros(string result)
{
    // removing trailing zero from sum, substract
    // eg 01 -> 1, 0001221022 -> 1221022, 00000-> 0
    int id = 0;
    for (id = 0; id < result.length(); id++)
    {
        if (result[id] != '0')
        {
            //cout << id << endl;
            break;
        }
    }
    //cout << "HERE " << result << endl;
    if (id == result.length())
        result = "0";
    else result = result.substr(id, result.length());
    return result;
}
string sum(string a, string b)
{
    // sum any big string "numbers"
    int cf = 0, temp=0;
    stack<int> st;
    string result = "";
    int size = max(a.length(), b.length());
    int lenA = a.length()-1;
    int lenB = b.length()-1;
    for (int i = 0; i < size;i++) // chua dung <-
    {
        if (lenA>= i && lenB >= i)
        {
            temp = a[lenA-i] -'0'+ b[lenB-i]-'0';
        }
        else if(lenA >= i)
        {
            temp = a[lenA - i] - '0';
        }
        else if (lenB>= i)
        {
            temp = b[lenB - i] - '0';
        }
        temp += cf;
        if (temp >= 10)
        {
            st.push(temp % 10);
            cf = temp / 10;
        }
        else
        {
            st.push(temp);
            cf = 0;
        }
        temp = 0;
    }
    if (cf != 0)
        st.push(cf);
    while (!st.empty())
    {
        result += to_string(st.top());
        st.pop();
    }
    result = checkTrailingZeros(result);
    return result;
}
string substract(string a, string b) {
    // substract any big string "numbers" 
    // only works for a value > b value
    int cf = 0;
    stack<int> st;
    string result = "";
    int size = max(a.length(), b.length());
    int lenA = a.length() - 1;
    int lenB = b.length() - 1;
    //cout << a << "---" << b << endl;
    for (int i = 0; i < size; i++)
    {
        int temp=0;
        if (lenB >= i)
            temp = b[lenB - i] - '0';
        temp += cf;
        if (temp > a[lenA - i] - '0')
        {
            st.push( (a[lenA - i] - '0') + 10 - temp);
            cf = 1;
        }
        else
        {
            st.push(a[lenA - i] - '0'-  temp);
            cf = 0;
        }
    }
    while (!st.empty())
    {
        result += to_string(st.top());
        st.pop();
    }
    //cout << result << endl;
    result = checkTrailingZeros(result);
    return result;
}

string decimal(int n, string x)
{
    // create 10^n*x
    for (int i = 0; i < n; ++i)
    {
        x += "0";
    }
    return x;
}
string quotient(string x, string p)
{
    // get x/10^n
    return x.substr(0, p.length() - 1);
}
string remainder(string x   , string p)
{
    // get x%10^n
    return x.substr(p.length()-1, x.length()-1);
}
string same_length(string a, int size)
{
    // adding trailing zeros to have equal length
    int d = size - a.length();
    if (d > 0)
    {
        return a.insert(0, d, '0');
    }
    else return a;
}
// karatsuba multiplication
string karatsuba(string x, string y)
{
    if (x.length() < 2 && y.length() < 2) // base case
    {
        long long res = stoll(x) * stoll(y);
        return to_string(res);
    }
    int size = fmax(x.length(), y.length());
    x=same_length(x, size);
    y=same_length(y, size);
    int n1 = (int)ceil(x.length() / 2.0);
    int n2 = (int)ceil(y.length() / 2.0);
    string p1 = decimal(n1,"1");//10n1
    string p2 = decimal(n2,"1");//10n2
    string a = quotient(x,p1);
    string b = remainder(x, p1);
    
    string c = quotient(y,p2);
    string d = remainder(y, p2);


    //// Recur until base case
    string ac = karatsuba(a, c);
    string bd = karatsuba(b, d);
    //(a+b)(c+d)-ac-bd
    string e = substract(substract(karatsuba(sum(a , b), sum(c , d)),ac),bd) ;

    //// if n is odd => size and n1 should minus 1
    if (size % 2 != 0)
    {

        size -= 1;
        n1 -= 1;
    }
    return sum(sum(decimal(size, ac), decimal(n1, e)), bd);
}
void multiply64digits(string a, string b)
{
    if (a.length() < 9 && b.length() < 9)
        cout << stoll(a) * stoll(b) << endl;
    else
        cout << karatsuba(a, b) << endl;
    
}