package Arrays;

import java.lang.reflect.Array;
import java.util.*;


public class Easy {
    public static void main(String[] args) {

    }

    //Find second largest element in an array
    /*
    Example:
array length = 3
input array={2,2,2}
return -1 as only (2) same value exist
if input array ={2,1,2}
return 1
     */
    //METHOD -1 ->SORTING   O(N LOG (N))   BRUTE

    public static int sec_lar_ele1(Integer[] arr) {
        int ele = -1;
        Arrays.sort(arr, Collections.reverseOrder());
        for (int i = 1; i < arr.length; i++) {
            if (!arr[i].equals(arr[0])) {
                return arr[i];
            }
        }
        return ele;
    }

    //    Meth 2 - Two Pass "    O(2 * N)    BETTER
    public static int sec_lar_ele(int[] arr) {
        if (arr == null || arr.length < 2) return -1;

        int max = Arrays.stream(arr).max().getAsInt();
        int sec_max = Integer.MIN_VALUE;

        for (int i = 0; i < arr.length; i++) {
            int ele = arr[i];
            if (ele < max && ele > sec_max) {
                sec_max = ele;
            }

        }
        return sec_max == Integer.MIN_VALUE ? -1 : sec_max;
    }

    //    Meth 3 - Single-Pass Linear Scan Approach"    O(N)   OPTIMAL
    public static int sec_lar_ele3(int[] arr) {
        if (arr == null || arr.length < 2) return -1;
        int max = Integer.MIN_VALUE;
        int sec_max = Integer.MIN_VALUE;
        for (int ele : arr) {
            if (ele > max) {
                max = ele;
            } else if (ele > sec_max && ele < max) {
                sec_max = ele;
            }
        }
/*        for second_smallest

             if (ele < min) {
                min = ele;
            } else if (ele < sec_min && ele > min) {
                sec_max = ele;
            }
  */
        return sec_max == Integer.MIN_VALUE ? -1 : sec_max;
    }

    //Check if the array is sorted
    public static boolean is_sorted(int[] arr) {
        int len = arr.length;
        for (int i = 0; i < len - 1; i++) {
            if (arr[i] > arr[i + 1]) return false;
        }
        return true;
    }
//    Remove duplicates from a sorted array returns count

    public static int rem_dup(int[] arr) {
// [1, 1, 2 , 2, 2, 3, 3]
        int i = 0, j = i + 1, len = arr.length;
        while (j < len) {
            if (arr[i] == arr[j]) {
                j++;
            } else if (arr[i] != arr[j]) {
                arr[i + 1] = arr[j];
                i++;
                j++;
            }
        }
//        fill underscores
        for (int k = i + 1; k < len; k++) {
            arr[k] = '_';
        }

        return i + 1;
    }
/* SHORT CODE
        int i = 0;
        //slide and place at  i
        for(int ele:arr){
         if(i<1 || ele>arr[i-1]){
             arr[i++] = ele;
         }
        }
    return i;

 */

    //    At most 2 unique are allowed
    public static int rem_dup_AtMost2(int[] arr) {
        int i = 0;
        for (int ele : arr) {
            if (i < 2 || ele > arr[i - 2]) {
                arr[i] = ele;
                i++;
            }

        }
        return i;
    }

    //    Remove duplicates from an unsorted array
    public static int rem_dup_unsorted(int[] arr) {


        return 0;
    }


    //    Rotate the array by k place
//    Brute force
    public void rotate_brute(int[] arr, int k) {
//  ALGO -->   shifting and store
/*        for k = 1
        int temp = arr[0],n= arr.length; //store

        for (int i = 1; i <n ; i++) { //shift
            arr[i-1] = arr[i];
        }
//        replace
        arr[n-1] = temp;
  */
        int n = arr.length;
        k = k % n;
        int[] temp = new int[k];
        for (int i = 0; i < k; i++) { //store
            temp[i] = arr[i];
        }
        for (int i = k; i < n; i++) { //shift
            arr[i - k] = arr[i];
        }
        //replace

        for (int i = n - k; i < n; i++) {
            arr[i] = temp[i - (n - k)];
        }
//T.C--> O(N + K) S.C -> O(K)
    }

    //    Optimized
//    T.C--> O(2N) SC-> O(1) OF EXTRA SPACE
    public void rotate(int[] arr, int k) {
        int n = arr.length;
        k = k % n;
        if (k == 0) return;//no rotation required
        k_rev(arr, 0, n - k - 1);
        k_rev(arr, n - k, n - 1);
        k_rev(arr, 0, n - 1);
    }

    //    Customized reverse function
    public static void k_rev(int[] arr, int i, int j) {
        while (i < j) {
            swap(arr, i, j);
            i++;
            j--;
        }

    }

    public static void swap(int[] arr, int i, int j) {
        int tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    //    MOVE ZEROES
    public static void move_zeroes(int[] arr) {
// ALGO-> shift and replace
//        [1 0 0 3 12]
        int i = 0;
//pick each element
        for (int ele : arr) {

            if (ele != 0) {
                arr[i++] = ele;
            }
        }
        for (int j = i; j < arr.length; j++) {
            arr[j] = 0;
        }

/*        code 2 --> striver code
//  locate first 0 and point j
  int j = -1;
  for(int i = 0;i<n;i++){
  if(arr[i]==0){
  j=i;
  break;
  }
  }
  if(j==-1)return ; //no zeroes found
// shifting
  for(int i = j+1;i<n;i++){
  if(arr[i]!=0){
        swap(arr[i],arr[j]);
        j++;
    }

*/
    }

    //    Union of two sorted arrays
    //Brute force --> Using set
    public static int[] Union(int[] arr1, int[] arr2) {
        TreeSet<Integer> set = new TreeSet<>();
        for (int ele : arr1) set.add(ele); // n log(n)
        for (int ele : arr2) set.add(ele);  // m log(m)
        int[] res = new int[set.size()];
        int i = 0;
        for (int ele : set) {    // n + m
            res[i++] = ele;
        }
        return res;

        //T.C --> NLOGN + MLOGN + (N+M)
//        S.C --> O(N+M) OF external  SET
//Union set space used to return answer specifically mention this to interviewer
    }

    //    Optimized --> USING TWO POINTERS
    public static int[] Union_optimized(int[] arr1, int[] arr2) {
        int n = arr1.length, m = arr2.length, i = 0, j = 0;
        List<Integer> ls = new ArrayList<>();
//        Move Pointers
        while (i < n && j < m) {
            if (arr1[i] < arr2[j]) {
                if (ls.isEmpty() || arr1[i] != ls.get(ls.size() - 1)) {
                    ls.add(arr1[i]);
                }
                i++;
            } else if (arr1[i] > arr2[j]) {
                if (ls.isEmpty() || arr2[j] != ls.get(ls.size() - 1)) {
                    ls.add(arr2[j]);
                }
                j++;
            } else { //arr1[i] = = arr2[j]
                if (ls.isEmpty() || arr1[i] != ls.get(ls.size() - 1)) {  //checking Uniqueness
                    ls.add(arr1[i]);
                }
                i++;
                j++;
            }
        }
//Add remaining elements
        while (i < n) {
            if (arr1[i] != ls.get(ls.size() - 1)) {
                ls.add(arr1[i]);
            }
            i++;
        }
        while (j < m) {
            if (arr2[j] != ls.get(ls.size() - 1)) {
                ls.add(arr2[j]);
            }
            j++;
        }


//        Preparing result array
        int[] res = new int[ls.size()];
        for (int l = 0; l < ls.size(); l++) {
            res[l] = ls.get(l);
        }
        return res;
//        T.C --> O(N + M) , SC --> O(N + M) FOR RETURNING THE ANS NOT USED IN MY ALGO
    }

    //INTERSECTION OF TWO SORTED ARRAYS
//    Brute force --> compare each element

    //    All the elements present in both repetition is allowed
    public static int[] InterSection(int[] arr1, int[] arr2) {
        int n = arr1.length, m = arr2.length, i = 0, j = 0;
        List<Integer> ls = new ArrayList<>();
//Move Pointers

        while (i < n && j < m) {
            if (arr1[i] != arr2[j]) {
                if (arr1[i] < arr2[j]) {
                    i++;
                } else {
                    j++;
                }
            } else {
                ls.add(arr1[i]);
                i++;
                j++;
            }
        }

        int[] res = new int[ls.size()];
        for (int k = 0; k < ls.size(); k++) {
            res[k] = ls.get(k);
        }
        return res;
    } //T.C --> O(N+M) , S.C ->O(1) or mention O(K)  WHERE K = MIN(N,M) space required only to return result

    //Missing Number
    public static int missing_number_Brute(int[] arr) {
        //T.C --> o(N*N)
        int n = arr.length;
        boolean flag = true;
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j < n; j++) {
                if (arr[j] == i) {
                    flag = false;
                    break;
                }
            }
            if (flag) return i;
        }
        return -1; //no element found
    }

    public static int missing_number_Better(int[] arr) {
        //  Hashing    T.C --> o(N)  , S.C --> O(N)
        int n = arr.length;
        int[] hash = new int[n + 1];
        Arrays.fill(hash, -1);
        for (int i = 0; i < n; i++) {
            hash[arr[i]] = i;
        }
        for (int i = 0; i < hash.length; i++) {
            if (hash[i] == -1) return i;
        }
        return -1;
    }

    public static int missing_number_Optimal1(int[] arr) {
//    Algo --> Summative   //i will overflow for the number like 10^5
        int sum = 0, n = arr.length;
        for (int ele : arr) sum += ele;
        return ((n * (n + 1) / 2) - sum);
    }

    public int missingNumber_optimal2(int[] arr) {
        int n = arr.length;
        int xor1 = 0, xor2 = 0;

        // XOR all elements in the array
        for (int i = 0; i < n; i++) {
            xor2 ^= arr[i];
        }

        // XOR all numbers from 0 to n
        for (int i = 0; i <= n; i++) {
            xor1 ^= i;
        }

        // The missing number is the result of XORing xor1 and xor2
        return xor1 ^ xor2;
    }
    public int findMaxConsecutiveOnes(int[]arr) {
        int sum = 0,max_ele = 0;
        for(int ele:arr){
            if(ele==1){
                sum++;
            }else{
                max_ele = Math.max(max_ele,sum);
                sum = 0;
            }
        }
        return       max_ele = Math.max(max_ele,sum);
    }

    //    Find the element that appear once
    public int singleNumber(int[] arr) {
/*Algo
Brute --> loop each element and count   T.C--> O(n^2)
Better --> Hashing   1. Array hashing  get max and hash array and check cnt
its fails when number is -ve or too big like 10^9 or big
2. Using Map data structure - T.C  --> N(LOG M) depends on ordered and Unordered map+ (n/2 +1)
[loop to search in map]
M is the size of the map which is   (n/2 +1).
S.C --> (N/2 +1)  used by the map

optimal --> using xor

 */

        return -1;
    }
}