#include <string.h>

template <typename T>
void colToRowOrder(T* A, int rows, int cols){
  T* newStorage = new T[rows * cols];

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      newStorage[i * cols + j] = A[j * rows + i];
    }
  }

  memcpy(A, newStorage, sizeof(T) * rows * cols);
  delete[] newStorage;
}

template <typename T>
void rowToColOrder(T* A, int rows, int cols){
  T* newStorage = new T[rows * cols];

  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      newStorage[j * rows + i] = A[i * cols + j];
    }
  }

  memcpy(A, newStorage, sizeof(T) * rows * cols);
  delete[] newStorage;
}
