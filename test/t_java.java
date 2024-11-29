import nlopt.*;

public class t_java {
  private static double myfunc(double[] x, double[] grad) {
    if (grad != null) {
      grad[0] = 0.0;
      grad[1] = 0.5 / Math.sqrt(x[1]);
    }
    return Math.sqrt(x[1]);
  }

  private static double myconstraint(double[] x, double[] grad, double a,
                                     double b) {
    if (grad != null) {
      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[1] = -1.0;
    }
    return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
  }

  public static void main(String[] args) {
    try {
      System.loadLibrary("nloptjni");
    } catch (UnsatisfiedLinkError e) {
      // Dependent libraries are not looked up in the java.library.path.
      // See: https://bugs.openjdk.org/browse/JDK-8213772
      // So they can end up not getting found. At least on Windows, this
      // workaround can make it work. (On GNU/Linux, you probably want to use
      // RPATH/RUNPATH or LD_LIBRARY_PATH instead.)
      // Note that this assumes a shared nlopt library. If nlopt is static, we
      // should not get here in the first place, nloptjni should load fine.
      System.loadLibrary("nlopt");
      System.loadLibrary("nloptjni");
    }
    Algorithm algo = args.length < 1 ? Algorithm.LD_MMA
      : Algorithm.swigToEnum(Integer.parseInt(args[0]));
    Opt opt = new Opt(algo, 2);
    System.out.println("algo: " + opt.getAlgorithmName());
    opt.setLowerBounds(new DoubleVector(Double.NEGATIVE_INFINITY, 1e-6));
    opt.setMinObjective(t_java::myfunc);
    opt.addInequalityConstraint((x, grad) -> myconstraint(x, grad, 2, 0), 1e-8);
    opt.addInequalityConstraint((x, grad) -> myconstraint(x, grad, -1, 1),
                                1e-8);
    opt.setXtolRel(1e-4);
    DoubleVector x0 = new DoubleVector(1.234, 5.678);
    DoubleVector x = opt.optimize(x0);
    double minf = opt.lastOptimumValue();
    System.out.println("optimum at " + x);
    System.out.println("minimum value: " + minf);
    System.out.println("result code: " + opt.lastOptimizeResult());
    System.out.println("nevals: " + opt.getNumevals());
    System.out.println("initial step: " + opt.getInitialStep(x0));
    assert Math.abs(minf - 0.544331) < 1e-3: "wrong optimum";
  }
}
