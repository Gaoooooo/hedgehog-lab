import { Mat } from './matrix';
import { lup, qr } from 'mathjs';

// Cholesky
export class Chol {
  L: Mat;
  constructor(A: Mat) {
    if (A.rows !== A.cols || A.rows === 0 || A.cols === 0) {
      throw new Error('Wrong dimension of matrix A.');
    }

    //dimension n
    const n = A.rows;

    //check for hermitian (note: does not account for complex/imaginary matrix values)
    // while this could be added to the later NaN check, this lets the user know what caused the NaN: positive definite, or hermitian errors
    if ( !(A.equals( A.T() ) ) ) {
      throw new Error('Non Hermitian; Cholesky decomposition does not apply.');
    }

    //matrix L
    let L = new Mat().zeros(n, n);


    //iteration
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < i + 1; k++) {
        let sum = 0;
        for (let j = 0; j < k; j++) {
          sum += L.val[i][j] * L.val[k][j];
        }

        if (i === k) {
          L.val[i][k] = Math.sqrt(A.val[i][i] - sum);
          //if giving NaN value (not positive definite), throw error
          if (isNaN(L.val[i][k])) {
            throw new Error('Not positive definite. Cholensky decomposition does not apply');
          }
        } else {
          L.val[i][k] = (1.0 / L.val[k][k]) * (A.val[i][k] - sum);
          if (isNaN(L.val[i][k])) {
            throw new Error('Not positive definite.');
          }
        }
      }
    }

    this.L = L;
  }
}

// LU
export class LU {
  L: Mat;
  U: Mat;
  P: Mat;
  constructor(A: Mat) {
    let result = lup(A.val);
    this.L = new Mat(result.L as any);
    this.U = new Mat(result.U as any);
    this.P = new Mat(result.P);
  }
}

// QR
export class QR {
  Q: Mat;
  R: Mat;
  constructor(A: Mat) {
    let result = qr(A.val);
    this.Q = new Mat(result.Q as any);
    this.R = new Mat(result.R as any);
  }
}
