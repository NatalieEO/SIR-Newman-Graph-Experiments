/* 
 * File:   utility.hpp
 * Author: maleq
 *
 * Created on November 14, 2008, 3:03 PM
 */

#ifndef _UTILITY_HPP
#define	_UTILITY_HPP

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <fstream>
#include <sys/time.h>
#include "Random.hpp"

using namespace std;

typedef  unsigned long int  ULI;


#define  Min(x,y) ((x)<(y) ? (x) : (y))
#define  Max(x,y) ((x)>(y) ? (x) : (y))
#define  Sqr(x) ((x)*(x))


/*-----------------------FindMax ------------------------------------------*/
template <typename T, typename SizeT>
T FindMax(T *d, SizeT  size) {
	T max = d[0];
	for (SizeT i=1; i<size; i++) 
		if (max < d[i])
			max = d[i];
	return max;		
}


/*------------------------ Round --------------------------------------------*/
inline int Round(double x)
{
	int ix = (int) x;
	if (x - ix < 0.5)
		return ix;
	return (ix+1);
}

/*--------------------- Swap the content of x and y ------------------------*/
template <typename T>
inline void Swap(T &x, T &y)
{
	T tmp = x;
	x = y;
	y = tmp;
}



/*---------------- Sort elements in d in ascending order -------------------*/
template <typename T, typename SizeT>
void Sort(T *d, SizeT  left, SizeT right)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if (d[i] < d[left]) {
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) Sort(d, left, k-1);
	if (k+1 < right) Sort(d, k+1, right);
	return;
}


/*---------------- Sort the records in d in ascending order ----------------*/
template <typename T, typename SizeT>
void RecSort(T *d, SizeT  left, SizeT right)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if (d[i].Key() < d[left].Key()) {
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) RecSort(d, left, k-1);
	if (k+1 < right) RecSort(d, k+1, right);
	return;
}



/*---------------- Sort elements in d in descending order ------------------*/
template <typename T, typename SizeT>
void SortD(T *d, SizeT  left, SizeT right)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if (d[i] > d[left]) {
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) SortD(d, left, k-1);
	if (k+1 < right) SortD(d, k+1, right);
	return;
}


/*----- Sort elements in d based on: first its label, then itself ----------*/
template <typename T, typename SizeT, typename LabelT>
void Sort(T *d, SizeT  left, SizeT right, LabelT *label)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if ( (label[d[i]] < label[d[left]]) 
			|| (label[d[i]]==label[d[left]] && d[i]<d[left]) )
		{
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) Sort(d, left, k-1, label);
	if (k+1 < right) Sort(d, k+1, right, label);
	return;
}



/*-------------- Sort elements in d based on label only --------------------*/
template <typename T, typename SizeT, typename LabelT>
void LSort(T *d, SizeT  left, SizeT right, LabelT *label)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if (label[d[i]] < label[d[left]]) {
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) LSort(d, left, k-1, label);
	if (k+1 < right) LSort(d, k+1, right, label);
	return;
}



/*-------------- Sort records in d based on label only --------------------*/
template <typename T, typename SizeT, typename LabelT>
void LRecSort(T *d, SizeT  left, SizeT right, LabelT *label)
{
	SizeT  i, k;
	
	for (k=left, i=left+1; i <= right; i++) {
		if (label[d[i].Key()] < label[d[left].Key()]) {
			k++;
			Swap(d[i], d[k]);
		}
	}
	Swap(d[k], d[left]);
	if (left+1 < k) LRecSort(d, left, k-1, label);
	if (k+1 < right) LRecSort(d, k+1, right, label);
	return;
}



/*--------------------- Print the element of an array ----------------------*/
template <typename ArrayT, typename T>
void PrintList(ArrayT *list, T size)
{
	cout << "(" << size << "): "; 
	for (T i=0; i<size; i++)
		cout << (int) list[i] << " ";
	cout << endl;
}


/*--------------------- Print an array of records ----------------------*/
template <typename ArrayT, typename T>
void PrintRecList(ArrayT *list, T size)
{
	cout << "(" << size << "): "; 
	for (T i=0; i<size; i++)
		list[i].Print();
}


/*--------------------- Print an array with marker i ----------------------*/
template <typename MarkerT, typename ArrayT, typename T>
void PrintList(MarkerT i, ArrayT *list, T size)
{
	cout << i << "(" << size << "): "; 
	for (T i=0; i<size; i++)
		cout << (int) list[i] << " ";
	cout << endl;
}

/*--------------------- Print an array with marker i ----------------------*/
template <typename MarkerT, typename ArrayT, typename T>
void PrintRecList(MarkerT i, ArrayT *list, T size)
{
	cout << i << " (" << size << "): "; 
	for (T i=0; i<size; i++)
		list[i].Print();
	cout << endl;
}


/*----------------------- Fill array a with value -------------------------*/
template <typename T, typename SizeT>
inline void Fill(T *a, SizeT size, T value)
{
	for (SizeT i=0; i<size; i++)
		a[i] = value;
}



/*---- MappedFill fills only the mapped elements (defined by map) of a ----*/
template <typename T, typename MapT, typename SizeT>
inline void MFill(T *a, MapT *map, SizeT size, T value)
{
	for (SizeT i=0; i<size; i++)
		a[map[i]] = value;
}



/*-------------------------- Copy array b to a -----------------------------*/
template <typename T, typename SizeT>
inline void Copy(T *a, T *b, SizeT size)
{
	for (SizeT i=0; i<size; i++)
		a[i] = b[i];
}



/*------- Find intersection of sorted arrays a and b, and store in a ------*/
template <typename T, typename SizeT>
inline SizeT Intersection(T *a, SizeT sa, T *b, SizeT sb)
{
	SizeT s, i, k;
	
	for (s=i=k=0; k<sb && i<sa; i++) {
		while (k<sb && b[k]<a[i]) k++;
		if (k<sb && b[k]==a[i]) a[s++] = a[i];
	}
	return s;
}



/*------- Find intersection of sorted arrays a and b, and store in c ------*/
template <typename T, typename SizeT>
inline SizeT Intersection(T *c, T *a, SizeT sa, T *b, SizeT sb)
{
	SizeT s, i, k;
	
	for (s=i=k=0; k<sb && i<sa; i++) {
		while (k<sb && b[k]<a[i]) k++;
		if (k<sb && b[k]==a[i]) c[s++] = a[i];
	}
	return s;
}


/*------------------ search item in UNSORTED array a ----------------------*/
template <typename T, typename SizeT>
inline bool LSearch(T *a, SizeT size, T item)
{
	for (SizeT i=0; i<size; i++) 
		if (a[i]==item) return true;
	return false;
}


/*------------- search item in UNSORTED array a of records ----------------*/
template <typename T, typename SizeT, typename KeyT>
inline bool LRecSearch(T *a, SizeT size, KeyT item)
{
	for (SizeT i=0; i<size; i++) 
		if (a[i].Key()==item) return true;
	return false;
}


/*-------- Binary search item in SORTED array a ---------------------------*/
template <typename T, typename SizeT>
inline SizeT BSearch(T *a, SizeT left, SizeT right, T item)
{
	SizeT mid;
	while (left <= right) {
		mid = (left + right) / 2;
		if (item == a[mid]) return (mid+1);
		if (item < a[mid]) right = mid - 1;
		else left = mid + 1;
	}
	return 0;
}



/*-------- Binary search item in array a SORTED in descending order -------*/
template <typename T, typename SizeT>
inline SizeT BSearchD(T *a, SizeT left, SizeT right, T item)
{
	SizeT mid;
	while (left <= right) {
		mid = (left + right) / 2;
		if (item == a[mid]) return (mid+1);
		if (item > a[mid]) right = mid - 1;
		else left = mid + 1;
	}
	return 0;
}



/*-------- Binary search item in SORTED array a ---------------------------*/
template <typename T, typename SizeT, typename KeyT>
inline bool BRecSearch(T *a, SizeT left, SizeT right, KeyT item)
{
	SizeT mid;
	while (left <= right) {
		mid = (left + right) / 2;
		if (item == a[mid].Key()) return true;
		if (item < a[mid].Key()) right = mid - 1;
		else left = mid + 1;
	}
	return false;
}



/*-------- Binary range search item in SORTED array a ---------------------*/
// returns 0 if item < a[0], returns i if a[i-1] <= item < a[i]
// input item should be less than a[right], returns right if item >= a[right]
template <typename T, typename SizeT>
inline SizeT BRSearch(T *a, SizeT left, SizeT right, T item)
{
	SizeT mid;
	while (left < right) {
		mid = (left + right) / 2;
		if (item < a[mid]) right = mid;
		else left = mid+1;
	}
	return left;
}


/*-- Binary range search LABELED item in SORTED (based on label) array a ---*/
// returns left if item < a[left], returns i if a[i-1] <= item < a[i]
// input item should be less than a[right], returns right if item >= a[right]
template <typename T, typename SizeT, typename LT>
inline SizeT BRSearch(T *a, SizeT left, SizeT right, LT item, LT *label)
{
	SizeT mid;
	while (left < right) {
		mid = (left + right) / 2;
		if (item < label[a[mid]]) right = mid;
		else left = mid+1;
	}
	return left;
}


/*-- Binary range search LABELED item in SORTED (based on label) array a ---*/
// returns left if item < a[left], returns i if a[i-1] <= item < a[i]
// input item should be less than a[right], returns right if item >= a[right]
template <typename T, typename SizeT, typename LT>
inline SizeT BRRecSearch(T *a, SizeT left, SizeT right, LT item, LT *label)
{
	SizeT mid;
	while (left < right) {
		mid = (left + right) / 2;
		if (item < label[a[mid].Key()]) right = mid;
		else left = mid+1;
	}
	return left;
}


/*---------------- Remove item from array d -------------------------------*/
template <typename T, typename SizeT>
inline SizeT  Remove(T *d, SizeT  left, SizeT right, T item)
{
	SizeT  i, k;
	
	for (k=right+1, i=right; i >= left; i--) 
		if (d[i]==item)  d[i] = d[--k];
	return k;
}


/*---------------- Remove item from tail of array d -----------------------*/
template <typename T, typename SizeT>
inline SizeT  RemoveTail(T *d, SizeT  left, SizeT right, T item)
{
	SizeT  i;
	
	for (i=right; i>left && d[i]==item; i--);
	if (d[i]==0)
		return i;
	return (i+1);
}


/*---------------- FreeMem frees memory pointed by mem --------------------*/
template <typename T>
inline void  FreeMem(T &mem)
{
	if (mem)
		delete [] mem;
	mem = NULL; 
}



/*------- FreeAll frees memory pointed by itself and its element -----------*/
template <typename T, typename SizeT>
inline void  FreeAll(T &mem, SizeT s)
{
	if (mem) {
		for (SizeT i=0; i<s; i++)
			delete [] mem[i];
		delete [] mem;
	}
	mem = NULL; 
}


/*--------------- Fact returns the factorial of n -------------------------*/
template <typename T>
inline T  Factorial(T n)
{
	T result=1, i;
	for (i=2; i<=n; i++)
		result *= i;
	return result;
}


/*--------------- Skip next white spaces in a file ------------------------*/

template <typename FilePtr>
inline void SkipWhiteSpaces(FilePtr &fp)
{
	char c = fp.peek();
	while (isspace(c)) {
		fp.ignore();
		c = fp.peek();
	}
}



/*-- Select a random subset dest of sn distinct samples from src of size n --*/
/* Returns the number of samples selected, which is min(sn, n) */
/* caller is responsible to call srand or similar function */
template <typename DataT>
inline int RandomSubset(DataT *dest, int sn, DataT *src, int n)
{    
    if (n <= sn) {
        Copy (dest, src, n);
        return n;
    }
    
    int idx;
    DataT *tsrc = new DataT[n];
    Copy(tsrc, src, n);
    
    for (int i=0; i<sn; i++) {
        idx = Uniform(0,sn-i-1);
        dest[i] = tsrc[idx];
        tsrc[idx] = tsrc[sn-i-1];
    }
    delete [] tsrc;
    return sn;
}


/*-- Select a random subset dest of sn distinct samples from [0,n-1] --*/
/* Returns the number of samples selected, which is min(sn, n) */
/* caller is responsible to call srand or similar function */
template <typename DataT>
inline int RandomSubset(DataT *dest, int sn, int n)
{    
    int *tsrc = new int[n];
    for (int i=0; i<n; i++) 
        tsrc[i] = i;

    if (n <= sn) {
        Copy (dest, tsrc, n);
        delete [] tsrc;
        return n;
    }
    
    int idx;    
    for (int i=0; i<sn; i++) {
        idx = Uniform(0, sn-i-1);
        dest[i] = tsrc[idx];
        tsrc[idx] = tsrc[sn-i-1];
    }
    delete [] tsrc;
    return sn;
}



/*void PrintGalibCredit()
{
	cout << endl;
	system("cat /vbi/projects/CINET/Bin/galib_credit.txt");
	cout << endl;
}
*/

void PrintGalibCredit(const char *path)
{
	char command[1024];
	strcpy(command, "cat ");
	
	int i, k, l;
	for (i=strlen(path)-1; i>=0 && path[i]!='/'; i--);
	l = strlen(command);
	for (k=0; k<=i; k++)
		command[l+k] = path[k];
	command[l+k] = 0;
	strcat(command, "galib_credit.txt");
	cout << endl;
	system(command);
	cout << endl;
}


/*------- Calculate the difference of two time in "micro-second" ----------*/

long TimeDifference(struct timeval &start, struct timeval &end)
{
	struct timeval diff;
	diff.tv_sec  = end.tv_sec  - start.tv_sec;
	diff.tv_usec = end.tv_usec - start.tv_usec;
	while (diff.tv_usec < 0)
	{
		diff.tv_usec += 1000000;
		diff.tv_sec -= 1;
	}
        
	return 1000000L * diff.tv_sec + diff.tv_usec;	// returning time difference in "micro-second"
}



#endif	/* _UTILITY_HPP */

