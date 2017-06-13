import java.util.*;

// ===================================================
// Vector Class - Vectors in Java
// ===================================================
public class Vector {
    private double X, Y, Z;

    // Constructors
    public Vector(double x, double y, double z) {
	X = x; Y = y; Z = z;
    }

    // Accessors and Mutators
    public double getMagnitude() {
	return Math.sqrt(X*X + Y*Y + Z*Z);
    }
    public double getX() {return X;}
    public double getY() {return Y;}
    public double getZ() {return Z;}
    public Vector getCopy() {return new Vector(X,Y,Z);}
    public Vector getNormalized() {
	Vector v = getCopy();
	return v.normalize();
    }

    public double setX(double x) {
	double tmp = X;	X = x; return tmp;
    }
    public double setY(double y) {
	double tmp = Y;	Y = y; return tmp;
    }
    public double setZ(double z) {
	double tmp = Z;	Z = z; return tmp;
    }

    // Functions
    public Vector normalize() {
	double m = getMagnitude();
	if (m != 0) {
	    X /= m; Y /= m; Z /= m;
	}
	return this;
    }
    public double dot(Vector V) {
	return X * V.getX() + Y * V.getY() + Z * V.getZ();
    }
    public Vector add(Vector V) {
	return new Vector(X + V.getX(),
			  Y + V.getY(),
			  Z + V.getZ());
    }
    public Vector subtract(Vector V) {
	return new Vector(X - V.getX(),
			  Y - V.getY(),
			  Z - V.getZ());
    }
    public Vector scale(double s) {
	return new Vector(X * s,
			  Y * s, 
			  Z * s);
    }
}
