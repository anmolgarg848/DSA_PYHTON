package Arrays;

import java.util.Arrays;

public class Basic {
    public static void main(String[] args) {
        int arr[] = new int[10];
        for (int i = 0; i < 10; i++) {
            arr[i] = i;
        }
        System.out.println(Arrays.toString(arr));
        rev(arr);
        System.out.println(Arrays.toString(arr));
        System.out.println(fib(6));
        cnt1("avcdaaazzczbbcd");
//        int[] ar = {2, 5, 1, -3, -8, -1};
        int[] ar = {1,2 ,3,4,5};
//        selection_sort(ar);
//        Bubble_sort(ar);
        Insertion_sort(ar);
        System.out.println(Arrays.toString(ar));
    }

    //REVERSE AN ARRAY
    public static void rev(int[] arr) {
//        Two pointer swapping
        int len = arr.length;
        int i = 0, j = len - 1;
        while (i < j) {
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }

    }

    //    T.C-> O(N)
    public static int fib(int n) {
        int i = 0;
        int j = 1;

        while (n > 0) {
            int ans = i + j;
            i = j;
            j = ans;
            n--;
        }
        return i;
    }

    //Count elements of an array
//    Number Hashing
    public static int cnt(int[] arr) {
        int len = arr.length;
        int maxValue = Arrays.stream(arr).max().getAsInt();
        int[] tmp = new int[maxValue + 1];
        int res = 0;
        for (int i = 0; i < len; i++) {
            tmp[arr[i]]++;
        }
        return res;
    }

    //    Count frequency of a String
//    Character Hashing
    public static void cnt1(String str) {
        int len = str.length();
        int[] arr = new int[26];
        for (int i = 0; i < len; i++) {
            arr[str.charAt(i) - 'a']++;
        }
        System.out.println(Arrays.toString(arr));
    }

    //    Sorting Algorithms
    public static void selection_sort(int[] arr) {
//ALGO-> SELECT MINIMUM ELEMENT AND SWAP
//        T.C --> N, N-1 , N-2 , ....1 = (N*(N+1))/2 = O(N"2)
        int len = arr.length;
        for (int i = 0; i < len - 1; i++) {
            int min_idx = get_min_idx(arr, i, len - 1); // Find the index of the minimum element
            swap(arr, min_idx, i); // Swap the found minimum element with the first element
        }
    }

    public static void swap(int[] arr, int i, int j) {
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    public static int get_min_idx(int[] arr, int i, int j) {
        int min_idx = i;
        for (int k = i + 1; k <= j; k++) {
            if (arr[k] < arr[min_idx]) { // Find the index of the minimum element
                min_idx = k;
            }
        }
        return min_idx;
    }

    public static void Bubble_sort(int[] arr) {
//ALGO--> Pushes the maximum to the last by adjacent swapping
//        T.C- O(N'2)
//BEST CASE using optimization O(n) boolean flag if no swap's happen in first go
        int len = arr.length;
        for (int i = 0; i < len-1; i++) {
            for (int j = 0; j <len-i-1 ; j++) {
                if (arr[j]>arr[j+1]){
                    swap(arr,j,j+1);
                }
            }
        }

    }
    public static void Insertion_sort(int[] arr) {
//ALGO --> TAKES AN ELEMENT AND PLACES IT IN ITS CORRECT POSITION
//
        int len = arr.length;
        for (int i = 0; i < len; i++) {
            for (int j = i; j >0&&arr[j-1]>arr[j] ; j--) {
                    swap(arr,j,j-1);
            }
        }
    }

}