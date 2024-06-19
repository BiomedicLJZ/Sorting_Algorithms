#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <queue>
#include <string>
#include <iterator>
#include <fstream>
#include "gtest/gtest.h"

// Function to generate random numbers
template <typename T>
std::vector<T> generateRandomNumbers(typename std::vector<T>::size_type n) {
    std::vector<T> arr(n);
    std::random_device rd;
    std::mt19937 gen(rd());

    if constexpr (std::is_integral<T>::value) {
        std::uniform_int_distribution<T> dis(1, 10000);
        for (typename std::vector<T>::size_type i = 0; i < n; i++) {
            arr[i] = dis(gen);
        }
    } else if constexpr (std::is_floating_point<T>::value) {
        std::uniform_real_distribution<T> dis(0.0, 1.0);
        for (typename std::vector<T>::size_type i = 0; i < n; i++) {
            arr[i] = dis(gen);
        }
    }

    return arr;
}

// Bubble Sort
void bubbleSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                std::swap(arr[j], arr[j + 1]);
            }
        }
    }
}

// Insertion Sort
void insertionSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

// Selection Sort
void selectionSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = 0; i < n - 1; i++) {
        int minIndex = i;
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[minIndex]) {
                minIndex = j;
            }
        }
        std::swap(arr[i], arr[minIndex]);
    }
}

// Shell Sort
void shellSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int gap = n / 2; gap > 0; gap /= 2) {
        for (int i = gap; i < n; i++) {
            int temp = arr[i];
            int j;
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap) {
                arr[j] = arr[j - gap];
            }
            arr[j] = temp;
        }
    }
}

// Heapify function used in Heap Sort
void heapify(std::vector<int>& arr, int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && arr[left] > arr[largest]) {
        largest = left;
    }
    if (right < n && arr[right] > arr[largest]) {
        largest = right;
    }
    if (largest != i) {
        std::swap(arr[i], arr[largest]);
        heapify(arr, n, largest);
    }
}

// Heap Sort
void heapSort(std::vector<int>& arr) {
    int n = arr.size();
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(arr, n, i);
    }
    for (int i = n - 1; i >= 0; i--) {
        std::swap(arr[0], arr[i]);
        heapify(arr, i, 0);
    }
}

// Merge function used in Merge Sort
void merge(std::vector<int>& arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;
    std::vector<int> L(n1), R(n2);
    for (int i = 0; i < n1; i++) {
        L[i] = arr[l + i];
    }
    for (int j = 0; j < n2; j++) {
        R[j] = arr[m + 1 + j];
    }
    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k++] = L[i++];
        } else {
            arr[k++] = R[j++];
        }
    }
    while (i < n1) {
        arr[k++] = L[i++];
    }
    while (j < n2) {
        arr[k++] = R[j++];
    }
}

// Merge Sort
void mergeSort(std::vector<int>& arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}

// Partition function used in Quick Sort
int partition(std::vector<int>& arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Quick Sort
void quickSort(std::vector<int>& arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

// Bucket Sort
void bucketSort(std::vector<float>& arr) {
    int n = arr.size();
    std::vector<std::vector<float>> buckets(n);
    for (int i = 0; i < n; i++) {
        int bi = n * arr[i];
        buckets[bi].push_back(arr[i]);
    }
    for (int i = 0; i < n; i++) {
        std::sort(buckets[i].begin(), buckets[i].end());
    }
    int index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < buckets[i].size(); j++) {
            arr[index++] = buckets[i][j];
        }
    }
}

// Counting Sort used in Radix Sort
void countingSort(std::vector<int>& arr, int exp) {
    int n = arr.size();
    std::vector<int> output(n);
    std::vector<int> count(10, 0);
    for (int i = 0; i < n; i++) {
        count[(arr[i] / exp) % 10]++;
    }
    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }
    for (int i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }
    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }
}


// Radix Sort
void radixSort(std::vector<int>& arr) {
    int max = *std::max_element(arr.begin(), arr.end());
    for (int exp = 1; max / exp > 0; exp *= 10) {
        countingSort(arr, exp);
    }
}

// Function to check if the array is sorted
template <typename T>
bool isSorted(const std::vector<T>& arr) {
    for (size_t i = 0; i < arr.size() - 1; i++) {
        if (arr[i] > arr[i + 1]) {
            return false;
        }
    }
    return true;
}
void splitAndSort(const std::string& inputFile, const std::string& tempFilePrefix, int chunkSize) {
    std::ifstream input(inputFile);
    int fileCount = 0;

    while (input) {
        std::vector<int> numbers;
        for (int i = 0; i < chunkSize && input; i++) {
            int number;
            input >> number;
            numbers.push_back(number);
        }

        std::sort(numbers.begin(), numbers.end());

        std::ofstream tempFile(tempFilePrefix + std::to_string(fileCount));
        for (int number : numbers) {
            tempFile << number << " ";
        }

        fileCount++;
    }
}

void merge(const std::string& outputFile, const std::string& tempFilePrefix, int fileCount) {
    std::ofstream output(outputFile);
    std::vector<std::ifstream> tempFiles(fileCount);

    for (int i = 0; i < fileCount; i++) {
        tempFiles[i].open(tempFilePrefix + std::to_string(i));
    }

    std::vector<int> currentNumbers(fileCount, 0);
    std::vector<bool> isFileEmpty(fileCount, false);
    int emptyFileCount = 0;

    for (int i = 0; i < fileCount; i++) {
        if (!(tempFiles[i] >> currentNumbers[i])) {
            isFileEmpty[i] = true;
            emptyFileCount++;
        }
    }

    while (emptyFileCount < fileCount) {
        int minIndex = -1;
        for (int i = 0; i < fileCount; i++) {
            if (!isFileEmpty[i] && (minIndex == -1 || currentNumbers[i] < currentNumbers[minIndex])) {
                minIndex = i;
            }
        }

        output << currentNumbers[minIndex] << " ";

        if (!(tempFiles[minIndex] >> currentNumbers[minIndex])) {
            isFileEmpty[minIndex] = true;
            emptyFileCount++;
        }
    }
}

void externalSort(const std::string& inputFile, const std::string& outputFile, int chunkSize) {
    std::string tempFilePrefix = "temp";
    splitAndSort(inputFile, tempFilePrefix, chunkSize);
    merge(outputFile, tempFilePrefix, chunkSize);
}


struct Element {
    int value;
    int fileIndex;

    Element(int value, int fileIndex) : value(value), fileIndex(fileIndex) {}
};

struct Compare {
    bool operator()(const Element& a, const Element& b) {
        return a.value > b.value;
    }
};

void externalMultiwayMerge(const std::string& outputFile, const std::string& tempFilePrefix, int fileCount) {
    std::priority_queue<Element, std::vector<Element>, Compare> minHeap;
    std::vector<std::ifstream> tempFiles(fileCount);
    std::vector<int> currentNumbers(fileCount, 0);
    std::vector<bool> isFileEmpty(fileCount, false);
    int emptyFileCount = 0;

    for (int i = 0; i < fileCount; i++) {
        tempFiles[i].open(tempFilePrefix + std::to_string(i));
        if (!(tempFiles[i] >> currentNumbers[i])) {
            isFileEmpty[i] = true;
            emptyFileCount++;
        } else {
            minHeap.push(Element(currentNumbers[i], i));
        }
    }

    std::ofstream output(outputFile);

    while (!minHeap.empty()) {
        Element minElement = minHeap.top();
        minHeap.pop();
        output << minElement.value << " ";

        if (!(tempFiles[minElement.fileIndex] >> currentNumbers[minElement.fileIndex])) {
            isFileEmpty[minElement.fileIndex] = true;
            emptyFileCount++;
        } else {
            minHeap.push(Element(currentNumbers[minElement.fileIndex], minElement.fileIndex));
        }
    }
}

void polyphaseMerge(const std::string& outputFile, const std::string& tempFilePrefix, int fileCount) {
    std::priority_queue<Element, std::vector<Element>, Compare> minHeap;
    std::vector<std::ifstream> tempFiles(fileCount);
    std::vector<int> currentNumbers(fileCount, 0);
    std::vector<bool> isFileEmpty(fileCount, false);
    int emptyFileCount = 0;

    for (int i = 0; i < fileCount; i++) {
        tempFiles[i].open(tempFilePrefix + std::to_string(i));
        if (!(tempFiles[i] >> currentNumbers[i])) {
            isFileEmpty[i] = true;
            emptyFileCount++;
        } else {
            minHeap.push(Element(currentNumbers[i], i));
        }
    }

    std::ofstream output(outputFile);

    while (!minHeap.empty()) {
        Element minElement = minHeap.top();
        minHeap.pop();
        output << minElement.value << " ";

        if (!(tempFiles[minElement.fileIndex] >> currentNumbers[minElement.fileIndex])) {
            isFileEmpty[minElement.fileIndex] = true;
            emptyFileCount++;
            if (emptyFileCount == 1) {
                tempFiles[minElement.fileIndex].close();
                tempFiles[minElement.fileIndex].open(tempFilePrefix + std::to_string(fileCount));
                if (tempFiles[minElement.fileIndex] >> currentNumbers[minElement.fileIndex]) {
                    isFileEmpty[minElement.fileIndex] = false;
                    emptyFileCount--;
                    minHeap.push(Element(currentNumbers[minElement.fileIndex], minElement.fileIndex));
                }
            }
        } else {
            minHeap.push(Element(currentNumbers[minElement.fileIndex], minElement.fileIndex));
        }
    }
}

// Gtest for bubble sort
TEST(SortingTest, BubbleSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    bubbleSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Bubble Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for insertion sort
TEST(SortingTest, InsertionSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    insertionSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Insertion Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for selection sort
TEST(SortingTest, SelectionSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    selectionSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Selection Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for shell sort
TEST(SortingTest, ShellSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    shellSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Shell Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for heap sort
TEST(SortingTest, HeapSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    heapSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Heap Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for merge sort
TEST(SortingTest, MergeSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    mergeSort(arr, 0, arr.size() - 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Merge Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for quick sort
TEST(SortingTest, QuickSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    quickSort(arr, 0, arr.size() - 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Quick Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for bucket sort
TEST(SortingTest, BucketSortTest) {
    std::vector<float> arr = generateRandomNumbers<float>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    bucketSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Bucket Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for radix sort
TEST(SortingTest, RadixSortTest) {
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    auto start = std::chrono::high_resolution_clock::now();
    radixSort(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Radix Sort Time: " << elapsed.count() << "s\n";
    ASSERT_TRUE(isSorted(arr));
}

// Gtest for external sort
TEST(SortingTest, ExternalSortTest) {
    std::string inputFile = "input.txt";
    std::string outputFile = "output.txt";
    int chunkSize = 1000;  // Adjust this value based on the amount of available memory

    // Generate some random numbers and write them to the input file
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    std::ofstream input(inputFile);
    for (int number : arr) {
        input << number << " ";
    }

    // Perform the external sort
    auto start = std::chrono::high_resolution_clock::now();
    externalSort(inputFile, outputFile, chunkSize);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "External Sort Time: " << elapsed.count() << "s\n";

    // Read the sorted numbers from the output file and check if they are sorted
    std::ifstream output(outputFile);
    std::vector<int> sortedNumbers;
    int number;
    while (output >> number) {
        sortedNumbers.push_back(number);
    }
    ASSERT_TRUE(isSorted(sortedNumbers));
}

// Gtest for external multiway merge sort
TEST(SortingTest, ExternalMultiwayMergeSortTest) {
    std::string inputFile = "input.txt";
    std::string outputFile = "output.txt";
    std::string tempFilePrefix = "temp";
    int chunkSize = 1000;  // Adjust this value based on the amount of available memory

    // Generate some random numbers and write them to the input file
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    std::ofstream input(inputFile);
    for (int number : arr) {
        input << number << " ";
    }

    // Perform the external sort
    splitAndSort(inputFile, tempFilePrefix, chunkSize);

    // Perform the external multiway merge
    auto start = std::chrono::high_resolution_clock::now();
    externalMultiwayMerge(outputFile, tempFilePrefix, chunkSize);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "External Multiway Merge Sort Time: " << elapsed.count() << "s\n";

    // Read the sorted numbers from the output file and check if they are sorted
    std::ifstream output(outputFile);
    std::vector<int> sortedNumbers;
    int number;
    while (output >> number) {
        sortedNumbers.push_back(number);
    }
    ASSERT_TRUE(isSorted(sortedNumbers));
}

// Gtest for polyphase merge sort
TEST(SortingTest, PolyphaseMergeSortTest) {
    std::string inputFile = "input.txt";
    std::string outputFile = "output.txt";
    std::string tempFilePrefix = "temp";
    int chunkSize = 1000;  // Adjust this value based on the amount of available memory

    // Generate some random numbers and write them to the input file
    std::vector<int> arr = generateRandomNumbers<int>(100000);
    std::ofstream input(inputFile);
    for (int number : arr) {
        input << number << " ";
    }

    // Perform the external sort
    splitAndSort(inputFile, tempFilePrefix, chunkSize);

    // Perform the polyphase merge
    auto start = std::chrono::high_resolution_clock::now();
    polyphaseMerge(outputFile, tempFilePrefix, chunkSize);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Polyphase Merge Sort Time: " << elapsed.count() << "s\n";

    // Read the sorted numbers from the output file and check if they are sorted
    std::ifstream output(outputFile);
    std::vector<int> sortedNumbers;
    int number;
    while (output >> number) {
        sortedNumbers.push_back(number);
    }
    ASSERT_TRUE(isSorted(sortedNumbers));
}

int main() {
    testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}
