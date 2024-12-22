package Stacks_Queue;

import java.util.*;

public class SQ {
    //Strivers stacks and queue playlist
    public static void main(String[] args) {
        Stack<Integer> st = new Stack<>();
        Queue<Integer> q = new LinkedList<>();
        LinkedList<Integer> ls = new LinkedList<>();


    }
/*IMPLEMENT

//the only disadvantage is we need an constant size for below implementations
1. STACK USING ARRAYS    FNS --> PUSH , TOP , POP ,SIZE
2. QUEUE USING ARRAYS    FNS --> PUSH , TOP , POP ,SIZE

//Dynamic size
1. STACK USING LINKEDLIST    FNS --> PUSH , TOP , POP ,SIZE
2. QUEUE USING LINKEDLIST    FNS --> PUSH , TOP , POP ,SIZE

//Interesting problem
1. Implement stack using queue --> put element back in queue in reverse order
push will take O(n)  else will take o(1)

2. Implement queues using stack //Approach 1
NOTE : IN ALL QUEUE IMPLEMENTATIONS WE USES TWO VARIABLES START,END bcz In queue we are pushing in the
end and taking out from start
use 2 stack empty 1st into 2nd then push then empty 2nd into 1st
push will take O(2n)  else will take o(1)

//Approach 2 --> push will take less time but top and pop will take more time
simply push in stack 1
for top and pop if stack 2 is not empty then top and pop else transfer s1--->s2 then top and pop 

 */

    //BALANCED PARANTHESIS
    public static boolean Bal_Par(String str) {
        Stack<Character> st = new Stack<>();
        for (Character s : str.toCharArray()) {
            if (s == '{' || s == '[' || s == '(') {
                st.push(s);
            } else {
                if (st.isEmpty()) return false;
                char ch = st.pop();
                if (s == '}' && ch != '{' || s == ']' && ch != '[' || s == ')' && ch != '(') {
                    return false;
                }
            }
        }
        return st.isEmpty();
    }

    //Basic Calculator
    //Prefix, Infix, and Postfix Conversion
    /*
    Operator -->    ^ , * , / ,  + , -
    Operand -->   A-Z , a-z , 0-9
    Priority -->   3. ^  , 2. * , / , 1. +,- , -1. others
    prefix --> *+pq-mn
    infix --> (p+q)*(m-n)   //operators are between the operand
    postfix --> pq+mn-*


     */
    public static String infixTopostfix(String str) {
        //T.C --> O(N) + o(n)
        //S.C --> O(N) + o(n)
        StringBuilder sb = new StringBuilder();
        Stack<Character> st = new Stack<>();
        for (Character ch : str.toCharArray()) {
            //operands
            if (ch >= 'A' && ch <= 'Z' || ch >= 'a' && ch <= 'z' || ch >= '0' && ch <= '9') {
                //simply add in ans
                sb.append(ch);
            } else if (ch == '(') { //brackets
                st.add(ch);
            } else if (ch == ')') {
                //keep popping out
                while (!st.isEmpty() && st.peek() != '(') {
                    sb.append(st.pop());
                }
                st.pop();  //removing opening bracket
            } else { //operator
                //operator with more priority will be evaluate first in ans
                while (!st.isEmpty() && priority(st.peek()) >= priority(ch)) {
                    sb.append(st.pop());
                }
                st.push(ch);
            }
        }
        while (!st.isEmpty()) {
            sb.append(st.pop());
        }
        return sb.toString();
    }

    //    private static int priority(Character ch) {
//        return (ch == '^') ? 3 : (ch == '*' || ch == '/') ? 2 : (ch == '+' || ch == '-') ? 1 : -1;
//    }
    private static int priority(Character ch) {
        switch (ch) {
            case '^':
                return 3;  // Highest priority for exponentiation
            case '*':
            case '/':
                return 2;  // Medium priority for multiplication and division
            case '+':
            case '-':
                return 1;  // Lowest priority for addition and subtraction
            default:
                return -1; // For unrecognized operators
        }
    }

    //Infix to prefix conversion
    /* steps
    1. reverse the infix (make brackets vice-versa)
    2.do postfix conversion  //consider priority as strictly greater
    priority(st.peek()) > priority(ch)
    3. reverse that ans
     */
    public static String infixToprefix(String str) {
        //step 1 : reverse the infix (make brackets vice-versa)
        StringBuilder sb = new StringBuilder();
        for (int i = str.length() - 1; i >= 0; i--) {
            char ch = str.charAt(i);
            if (ch == '(') {
                sb.append(')');
            } else if (ch == ')') {
                sb.append('(');
            } else {
                sb.append(ch);
            }
        }
        //step 2 --> priority
        Stack<Character> st = new Stack<>();
        for (char ch : str.toCharArray()) {
            //operands
            if (ch >= 'A' && ch <= 'Z' || ch >= 'a' && ch <= 'z' || ch >= '0' && ch <= '9') {
                //simply add in ans
                sb.append(ch);
            } else if (ch == '(') { //brackets
                st.add(ch);
            } else if (ch == ')') {
                //keep popping out
                while (!st.isEmpty() && st.peek() != '(') {
                    sb.append(st.pop());
                }
                st.pop();  //removing opening bracket
            } else { //operator
                if (ch == '^') { //all operators with lower and equal priority shoud be added
                    //operator with more priority will be evaluate first in ans
                    while (!st.isEmpty() && priority(st.peek()) >= priority(ch)) {
                        sb.append(st.pop());
                    }
                } else { //other operators
                    while (!st.isEmpty() && priority(st.peek()) > priority(ch)) {
                        sb.append(st.pop());
                    }
                }
                st.push(ch);
            }
        }
        while (!st.isEmpty()) {
            sb.append(st.pop());
        }

//step 3 reverse again
        //t.c --> O(n/2)+O(n/2)+O(2n) ==> O(n)
        //s.c --> O(n) + O(n)
        return sb.reverse().toString();

    }

    /*POSTFIX TO INFIX CONVERSION
   example --> "AB-DE+F/*"
steps --> take last 2 elements from stack whenever encounters operator and place in between them and
wrap it up
     */
    public static String postfixToinfix(String str) {
        Stack<String> st = new Stack<>();
        for (int i = 0; i < str.length(); i++) {
            char ch = str.charAt(i);
            //operands
            if (ch >= 'A' && ch <= 'Z' || ch >= 'a' && ch <= 'z' || ch >= '0' && ch <= '9') {
                //simply add in stack
                st.add(ch + "");
            } else {//brackets not available in postfix
                //fetch last 2 operands
                String s2 = st.pop();
                String s1 = st.pop();
                st.add("(" + s1 + ch + s2 + ")");
            }
        }
        return st.pop();
    }

    //prefix to infix
    public static String prefixToinfix(String str) {
/*
"*+pq-mn"
iterating from the back
 */
        Stack<String> st = new Stack<>();
        for (int i = str.length(); i >= 0; i--) {
            char ch = str.charAt(i);
            //operands
            if (ch >= 'A' && ch <= 'Z' || ch >= 'a' && ch <= 'z' || ch >= '0' && ch <= '9') {
                //simply add in stack
                st.add(ch + "");
            } else {//brackets not available in postfix
                //fetch last 2 operands
                String s2 = st.pop();
                String s1 = st.pop();
                st.add("(" + s2 + ch + s1 + ")");
            }
        }
        return st.pop();
    }


    //postfix to prefix
    public static String postfixToprefix(String str) {
    /*
"AB-DE+F*'/'"
     */
        Stack<String> st = new Stack<>();
        for (int i = 0; i < str.length(); i++) {
            char ch = str.charAt(i);
            //operands
            if (ch >= 'A' && ch <= 'Z' || ch >= 'a' && ch <= 'z' || ch >= '0' && ch <= '9') {
                //simply add in stack
                st.add(ch + "");
            } else {//brackets not available in postfix
                //fetch last 2 operands
                String s2 = st.pop();
                String s1 = st.pop();
                st.add(ch + s1 + s2);
            }
        }
        return st.pop();

    }

    public static String prefixTOpostfix(String str) {
        //just iterate from back
        return "";
    }

    //Design minStack
    /*
    Approach 1--> store 2 elements in stack . S.C --> O(2N)
    Use pair class and keep track of min in top.second and return

    Optimized Approach
    Getting problem to backtrack to last minimum without storing it
    we solves this using mathematical calculation

    whenever we push value which is going to update my min then we insert modified value and replace
    minimum to the original value.
    2*val - prevMin =  newVal
//Take care
when you pop modified value then return min else return top
modified value will always less than min

     */
    class minStack {
        int min = Integer.MAX_VALUE, top = -1;
        Stack<Integer> st;

        minStack() {
            st = new Stack<>();
        }

        //push
        void push(int val) {
            if (st.isEmpty()) {
                min = val;
                st.add(val);
            } else {
                if (val < min) {
                    //Insert modified value
                    st.add(2 * val - min);
                    min = val;
                } else {
                    st.add(val);
                }
            }
        }

        //peak
        int peak() {
            if (st.isEmpty()) {
                System.out.println("error");
                return -1;
            } else {
                //check if value is modified or not
                if (st.peek() > min) { //not modified
                    return st.peek();
                } else {
                    return min;
                }
            }
        }

        void pop() {
            if (st.isEmpty()) {
                System.out.println("error");
            } else {
                int ele = st.pop();
                //check if value is modified or not
                if (ele < min) { // modified
                    //backtrack to prevMin
                    min = 2 * min - ele;
                }
            }
        }
    }

    /*MONOTONIC STACK --> stack that maintains its order in specific order

     */
    //NEXT GREATER ELEMENT
    public static int[] NGE(int[] arr) {
        //Intuition --> watching light poles as a human being
        //Stack will be in monotonic decreasing order
        int n = arr.length;
        int[] nge = new int[n];
        Stack<Integer> st = new Stack<>();
        //T.C --> O(2N)
        //Throughout the for loop at worst  stack will run for n times
        for (int i = n - 1; i >= 0; i--) {
            int ele = arr[i];
            while (!st.isEmpty() && st.peek() <= ele) st.pop();
            if (st.isEmpty()) {
                nge[i] = -1;
            } else {
                nge[i] = st.peek();
            }
            st.push(ele);
        }
        return nge;
    }

    public static int[] NGE_II(int[] arr) {
        //Intuition --> watching light poles as a human being
        //Stack will be in monotonic decreasing order
        int n = arr.length;
        int[] nge = new int[n];
        Stack<Integer> st = new Stack<>();
        //T.C --> O(2N)
        //Throughout the for loop at worst  stack will run for n times
        for (int i = 2 * n - 1; i >= 0; i--) {
            int ele = arr[i % n];
            while (!st.isEmpty() && st.peek() <= ele) st.pop();
            if (i < n) {
                if (st.isEmpty()) {
                    nge[i] = -1;
                } else {
                    nge[i] = st.peek();
                }
            }
            st.push(ele);
        }
        return nge;
    }
// Next smaller element from the left --same approach as above


    /*    Trapping rain water
    Intuition --> water logged on terris of each building so
    totalWaterLogged = min(LeftMax,RightMax)-heightOfCurBuilding

    Approach 1 -->preMax and sufMax array , T.C--> O(3*N),SC--> O(2N)
    Approach 2-> keep remeber the leftmax variable instead of array S.C-->(ON)
    Approach 3--> 2 way traversal and work on smaller building first and both pointers will always end
    up at largest building
     */
    public int trap(int[] arr) {
        int tot = 0, n = arr.length, l = 0, r = n - 1, lMax = 0, rMax = 0;
        while (l < r) {
            if (arr[l] <= arr[r]) { //left building
                if (lMax > arr[l]) tot += lMax - arr[l];
                else lMax = arr[l];
                l++;
            } else {//right building
                if (rMax > arr[r]) tot += rMax - arr[r];
                else rMax = arr[r];
                r--;
            }

        }
        return tot;
    }

    //Sum of SubArray Minimums --> T.C--> O(N^2)
    public int sumSubarrayMinimum(int[] arr) {
        int ans = 0, mod = (int) 1e9 + 7;
        for (int i = 0; i < arr.length; i++) {
            int min = arr[i];
            for (int j = i; j < arr.length; j++) {
                min = Math.min(min, arr[j]);
                ans += min % mod;
            }
        }
        return ans;
    }

    /*Optimized
    This gives the number of subarrays where arr[i] is the maximum, then we
    multiply it by arr[i] and add it to the total sum.

    The formula for counting subarrays including an element at index: (i+1) * (n-i)
    Intuition--> Individual contribution
    //edgecase --> consider whole array only once either from front or back [1,1]
     */
    public int sumSubarrayMinimum_optimized(int[] arr) {
        int n = arr.length, mod = (int) (1e9 + 7), tot = 0;
        int[] nseIdx = findNse(arr);  //next smaller elements
        int[] pseIdx = findPse(arr);  //previous smaller elements
        for (int i = 0; i < n; i++) {
            tot = (tot + (((nseIdx[i] - i) * (i - pseIdx[i])) * arr[i])) % mod;
        }
        return tot;
    }

    private int[] findNse(int[] arr) {
        int n = arr.length;
        int[] nxtIdx = new int[n];
        Stack<Integer> st = new Stack<>();
        for (int i = n - 1; i >= 0; i--) {
            int ele = arr[i];
            while (!st.isEmpty() && arr[st.peek()] >= ele) st.pop();
            if (st.isEmpty()) {
                nxtIdx[i] = n;
            } else {
                nxtIdx[i] = st.peek();
            }
            st.push(i);
        }
        return nxtIdx;

    }

    public static int[] findPse(int[] arr) {
        int n = arr.length;
        int[] prevIdx = new int[n];
        Stack<Integer> st = new Stack<>();
        for (int i = 0; i < n; i++) {
            int ele = arr[i];
            while (!st.isEmpty() && arr[st.peek()] > ele) st.pop(); //edge case handled
            if (st.isEmpty()) {
                prevIdx[i] = -1;
            } else {
                prevIdx[i] = st.peek();
            }
            st.push(i);
        }
        return prevIdx;
    }
    /*sum of subarray ranges
sum = tot_max-tot_min

     */

    public List<Integer> asteroidCollision(int[] arr) {
        List<Integer> ls = new ArrayList<>();
        for (int ele : arr) {
            if (ele > 0 || ls.isEmpty()) {  // Positive or stack is empty
                ls.add(ele);
            } else {  // Negative asteroid
                while (!ls.isEmpty() && ls.get(ls.size() - 1) > 0 && ls.get(ls.size() - 1) < -ele) {
                    ls.remove(ls.size() - 1);  // Destroy the smaller positive asteroid
                }
                if (!ls.isEmpty() && ls.get(ls.size() - 1) == -ele) {
                    ls.remove(ls.size() - 1);  // Both asteroids destroy each other
                } else if (ls.isEmpty() || ls.get(ls.size() - 1) < 0) {
                    ls.add(ele);  // No positive asteroids left to collide
                }
            }
        }
        return ls;
    }

    //Brute force --> Expand both sides and check
    public int largestRectangleArea(int[] heights) {
        //T.C --> O(N^2)
        int n = heights.length;
        int maxArea = heights[0];
        for (int i = 0; i < n; i++) {
            int right = i;
            //go right
            while (right < n && heights[i] <= heights[right]) right++;
            //go left
            int left = i;
            while (left >= 0 && heights[i] <= heights[left]) left--;
            int width = (right - i + 1) + (i - left);
            maxArea = Math.max(maxArea, width);
        }
        return maxArea;
        //optimized ---> use pse and nse approach
    }

    //remove k digit --> use stack
    //Maximal rectangle  --> refer histogram problem
//stock  span --> just use monotonic decreasing stack

    //Sliding window maximum
    /*
    Brute force
    T.C --> O((N-K)*K)
    optimized
    T.C--> O(2N) , S.C --> O(K)
     */
    public static int[] slidingWinMax(int[] arr, int k) {
        int n = arr.length;
        int[] max_arr = new int[n - k + 1];
        ArrayDeque<Integer> dq = new ArrayDeque<>();
        for (int i = 0; i < n; i++) {
            //maintain window of K elements
            if (!dq.isEmpty() && dq.peekFirst() <= (i - k)) {
                dq.pop();
            }
            while (!dq.isEmpty() && arr[i] >= arr[dq.peekLast()]) dq.pollLast();

            //add maximum element
            dq.addLast(i);
            //when first window starts add ans
            if (i >= (k - 1)) {
                max_arr[i - k + 1] = arr[dq.peekFirst()];
            }

        }
        return max_arr;
    }
//The Celebrity Problem
/*Brute Force
Two arrays know me and i know  approach
//Optimized --> two pointer top and bottom confirm and eliminate and then confirm
T.C --> O(2N)

 */



}
