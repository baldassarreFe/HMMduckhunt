public class HiddenMarkovModel {

	/** baum welch iterations */
	private static final long MAX_ITERATIONS = 50;

	public static double[] nextEmissions(double[][] A, double[][] B, double[] Pi) {
		return currentEmissions(B, nextStates(A, Pi));
	}

	public static double[] nextStates(double[][] A, double[] Pi) {
		return MatrixHelper.multiply(new double[][] { Pi }, A)[0];
	}

	public static double[] currentEmissions(double[][] B, double[] Pi) {
		int M = B[0].length;
		double[] result = new double[M];
		for (int currentEmission = 0; currentEmission < M; currentEmission++) {
			result[currentEmission] = MatrixHelper.multiply(Pi, MatrixHelper.getColumn(B, currentEmission));
		}
		return result;
	}

	public static double nextMostLikelyEmission(double[][] A, double[][] B, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(currentEmissions(B, nextStates(A, Pi)));
	}

	public static double nextMostLikelyState(double[][] A, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(MatrixHelper.multiply(new double[][] { Pi }, A)[0]);
	}

	public static double currentMostLikelyEmission(double[][] B, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(currentEmissions(B, Pi));
	}

	public static int[] estimateStateSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int T = O.length;
		int result[] = new int[T];
		@SuppressWarnings("unused")
		double probability = viterbi(O, A, B, Pi, result);
		return result;
	}

	public static double viterbi(int[] O, double[][] A, double[][] B, double[] Pi, int[] X) {
		int N = A.length;
		int T = O.length;

		double delta[][] = new double[T][N];
		int indexes[][] = new int[T][N];

		for (int i = 0; i < N; i++) {
			delta[0][i] = Pi[i] * B[i][O[0]];
			// useless because indexes[0] is never used
			indexes[0][i] = 0;
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				double max = 0;
				for (int j = 0; j < N; j++) {
					double temp = delta[t - 1][j] * A[j][i];
					if (temp > max) {
						max = temp;
						indexes[t][i] = j;
					}
				}
				delta[t][i] = max * B[i][O[t]];
			}
		}

		double probability = 0.0;
		int maxIndex = 0;
		for (int i = 0; i < N; i++) {
			if (probability < delta[T - 1][i]) {
				probability = delta[T - 1][i];
				maxIndex = i;
			}
		}

		for (int t = T - 1; t >= 0; t--) {
			X[t] = maxIndex;
			maxIndex = indexes[t][maxIndex];
		}

		return probability;
	}

	public static double viterbiNoLogs(int[] O, double[][] A, double[][] B, double[] Pi, int[] X) {
		int N = A.length;
		int T = O.length;

		double delta[][] = new double[T][N];
		int indexes[][] = new int[T][N];

		for (int i = 0; i < N; i++) {
			delta[0][i] = Math.log(Pi[i] * B[i][O[0]]);
			// useless because indexes[0] is never used
			indexes[0][i] = 0;
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				double max = Double.NEGATIVE_INFINITY;
				for (int j = 0; j < N; j++) {
					double temp = delta[t - 1][j] + Math.log(A[j][i]);
					if (temp > max) {
						max = temp;
						indexes[t][i] = j;
					}
				}
				delta[t][i] = max + Math.log(B[i][O[t]]);
			}
		}

		double probability = Double.NEGATIVE_INFINITY;
		int maxIndex = 0;
		for (int i = 0; i < N; i++) {
			if (probability < delta[T - 1][i]) {
				probability = delta[T - 1][i];
				maxIndex = i;
			}
		}

		for (int t = T - 1; t >= 0; t--) {
			X[t] = maxIndex;
			maxIndex = indexes[t][maxIndex];
		}

		return probability;
	}

	/**
	 * @param O
	 * @param A
	 * @param B
	 * @param Pi
	 * @return
	 */
	public static double[] baumWelch(int[] O, double[][] A, double[][] B, double[] Pi) {
		// initialization
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;
		int T = O.length;

		double oldProbLog;
		double logProb = Double.NEGATIVE_INFINITY;
		long iters = 0;
		do {
			oldProbLog = logProb;
			// alpha pass also gives log( P(O|lambda) )
			double[][] alpha = new double[T][N];
			double[] c = new double[T];
			logProb = alphaPass(O, A, B, Pi, alpha, c);

			double[][] beta = new double[T][N];
			betaPass(O, A, B, c, beta);

			double[][][] digamma = new double[T][N][N];
			double[][] gamma = new double[T][N];
			digammaGamma(O, A, B, alpha, beta, digamma, gamma);

			reEstimate(O, A, B, Pi, digamma, gamma);
		} while (++iters < MAX_ITERATIONS && logProb > oldProbLog);

		return new double[] { iters, logProb };
	}

	public static double alphaPass(int[] O, double[][] A, double[][] B, double[] Pi, double[][] alpha, double[] c) {
		int N = A.length;
		int T = O.length;

		c[0] = 0;
		for (int i = 0; i < N; i++) {
			alpha[0][i] = Pi[i] * B[i][O[0]];
			c[0] += alpha[0][i];
		}

		c[0] = 1.0 / c[0];
		for (int i = 0; i < N; i++) {
			alpha[0][i] *= c[0];
		}

		for (int t = 1; t < T; t++) {
			c[t] = 0;
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * A[j][i];
				}
				alpha[t][i] *= B[i][O[t]];
				c[t] += alpha[t][i];
			}
			c[t] = 1.0 / c[t];
			for (int i = 0; i < N; i++) {
				alpha[t][i] *= c[t];
			}
		}

		double logProb = 0;
		for (int t = 0; t < T; t++) {
			logProb += Math.log(c[t]);
		}
		return -logProb;
	}

	private static void alphaPassNoNormalization(int[] O, double[][] A, double[][] B, double[] Pi, double[][] alpha) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			alpha[0][i] = Pi[i] * B[i][O[0]];
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * A[j][i];
				}
				alpha[t][i] *= B[i][O[t]];
			}
		}
	}

	private static void betaPass(int[] O, double[][] A, double[][] B, double[] c, double[][] beta) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			beta[T - 1][i] = c[T - 1];
		}

		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
				beta[t][i] *= c[t];
			}
		}
	}

	private static void betaPassNoNormalization(int[] O, double[][] A, double[][] B, double[][] beta) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			beta[T - 1][i] = 1;
		}

		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
			}
		}
	}

	private static void digammaGamma(int[] O, double[][] A, double[][] B, double[][] alpha, double[][] beta,
			double[][][] digamma, double[][] gamma) {

		int N = A.length;
		int T = O.length;

		for (int t = 0; t < T - 1; t++) {
			double denom = 0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					denom += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
			}
			for (int i = 0; i < N; i++) {
				gamma[t][i] = 0;
				for (int j = 0; j < N; j++) {
					digamma[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j] / denom;
					gamma[t][i] += digamma[t][i][j];
				}
			}
		}

		double denom = 0;
		for (int i = 0; i < N; i++) {
			denom += alpha[T - 1][i];
		}
		for (int i = 0; i < N; i++) {
			gamma[T - 1][i] = alpha[T - 1][i] / denom;
		}
	}

	private static void digammaGammaNoNormalization(int[] O, double[][] A, double[][] B, double[][] alpha,
			double[][] beta, double[][][] digamma, double[][] gamma) {

		int N = A.length;
		int T = O.length;

		double denom = 0;
		for (int k = 0; k < N; k++) {
			denom += alpha[T - 1][k];
		}
		for (int t = 0; t < T - 1; t++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					digamma[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j] / denom;
				}
			}
		}
		for (int t = 0; t < T - 1; t++) {
			for (int i = 0; i < N; i++) {
				gamma[t][i] = 0;
				for (int j = 0; j < N; j++) {
					gamma[t][i] += digamma[t][i][j];
				}
			}
		}
	}

	private static void reEstimate(int[] O, double[][] A, double[][] B, double[] Pi, double[][][] digamma,
			double[][] gamma) {

		int N = A.length;
		int M = B[0].length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			// Pi
			Pi[i] = gamma[0][i];
		}
		for (int i = 0; i < N; i++) {
			// A
			for (int j = 0; j < N; j++) {
				double numerA = 0;
				double denomA = 0;
				for (int t = 0; t < T - 1; t++) {
					numerA += digamma[t][i][j];
					denomA += gamma[t][i];
				}
				A[i][j] = numerA / denomA;
			}
		}

		for (int i = 0; i < N; i++) {
			// B
			for (int j = 0; j < M; j++) {
				double numerB = 0;
				double denomB = 0;
				for (int t = 0; t < T; t++) {
					if (O[t] == j) {
						numerB += gamma[t][i];
					}
					denomB += gamma[t][i];
				}
				B[i][j] = numerB / denomB;
			}
		}
	}

	public static double probabilityOfEmissionSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int N = A.length;
		int T = O.length;
		double[][] alpha = new double[T][N];
		double[] c = new double[T];
		double result;

		alphaPass(O, A, B, Pi, alpha, c);
		result = 1;
		for (int t = 0; t < T; t++) {
			result *= c[t];
		}
		result = 1.0 / result;

		// alphaPassNoNormalization(O, A, B, Pi, alpha, c);
		// result = 0;
		// for (int i = 0; i < N; i++) {
		// result += alpha[T - 1][i];
		// }
		return result;
	}

	public static double probabilityLogOfEmissionSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int N = A.length;
		int T = O.length;
		double[][] alpha = new double[T][N];
		double[] c = new double[T];
		return alphaPass(O, A, B, Pi, alpha, c);
	}

}
