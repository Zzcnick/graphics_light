import java.io.*;
import java.util.*;

// ===================================================
// Pixel Class - Pixel Information
// ===================================================
public class Pixel {
    private int R;
    private int G;
    private int B;
    private double Z = Integer.MIN_VALUE;

    // Constructors
    public Pixel() { // White Pixel
	R = 255;
	G = 255;
	B = 255;
    }
    public Pixel(int R, int G, int B) {
	this.R = R;
	this.G = G;
	this.B = B;
    }
    public Pixel(int R, int G, int B, double Z) {
	this(R,G,B);
	this.Z = Z;
    }
    public Pixel(int[] RGB) {
	this(RGB[0], RGB[1], RGB[2]);
    }
    public Pixel(int[] RGB, double Z) {
	this(RGB[0], RGB[1], RGB[2]);
	this.Z = Z;
    }
    public Pixel(Pixel p) {
	int[] rgb = p.getRGB();
	R = rgb[0];
	G = rgb[1];
	B = rgb[2];
	Z = p.getZ();
    }
    public Pixel(int flag) {
	if (flag == 1) { // Random
	    R = (int)(Math.random() * 256);
	    G = (int)(Math.random() * 256);
	    B = (int)(Math.random() * 256);
	}
	if (flag == 2) { // Red Shades
	    R = (int)(Math.random() * 200 + 56);
	    G = (int)(R * 0.7);
	    B = (int)(R * 0.7);
	}
    }
    
    // Accessors and Mutators
    public int[] getRGB() {
	return new int[]{R, G, B};
    }
    public double getZ() {
	return Z;
    }
    
    public double setZ(double Z) {
	return this.Z = Z;
    }

    // Methods
    public int[] adjust(int dR, int dG, int dB) {
	int[] old = getRGB();
	R += dR; G += dG; B += dB;
	return old;
    }
    public boolean compareZ(Pixel p) {
	return Z >= p.getZ();
    }
    public Pixel copy() {
	return new Pixel(R, G, B);
    }

    // ToString Utility
    public String toString() {
	return "" + R + " " + G + " " + B + "\n";
    }
}
