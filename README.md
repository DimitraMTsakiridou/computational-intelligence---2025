# Computational Intelligence - 8th Semester

This repository contains four assignments completed for the **Computational Intelligence** course during my 8th semester.

## 📚 Assignments Overview

### 1. **Car Control using Fuzzy Logic Controller (FLC)**
- **Objective**: Design a fuzzy logic controller for vehicle obstacle avoidance
- **Key Components**:
  - Fuzzy sets for vertical distance (`dV`), horizontal distance (`dH`), and direction (`Θ`)
  - Mamdani inference system with max-min composition
  - Center of Area (COA) defuzzification
  - MATLAB implementation with FIS Editor
- **Features**: Vehicle navigation to target position (10, 3.2) with obstacle avoidance

### 2. **Classification using TSK Fuzzy Models**
- **Objective**: Solve classification problems using TSK fuzzy models
- **Datasets**: 
  - Haberman's Survival Dataset (306 instances, 3 features)
  - Epileptic Seizure Recognition Dataset (11,500 instances, 179 features)
- **Methodology**:
  - Subtractive Clustering for rule generation
  - Hybrid learning (backpropagation + least squares)
  - Feature selection and cross-validation
  - Performance metrics: Error matrix, OA, PA, UA, Kappa

### 3. **Regression using TSK Fuzzy Models**
- **Objective**: Model multivariate nonlinear functions using TSK models
- **Datasets**:
  - Airfoil Self-Noise Dataset (1,503 instances, 6 features)
  - Superconductivity Dataset (21,263 instances, 81 features)
- **Approach**:
  - Singleton and polynomial output functions
  - Bell-shaped membership functions
  - Grid search with 5-fold cross-validation
  - Performance metrics: RMSE, NMSE, NDEI, R²

### 4. **Satellite Attitude Control using Fuzzy Logic**
- **Objective**: Design fuzzy controllers for satellite orientation angle
- **Components**:
  - Fuzzy PI controller (FZ-PI)
  - Nine linguistic variables for error, error change, and control output
  - Mamdani inference with Center of Sums defuzzification
  - Comparison with linear PI controller
- **Features**: Reference tracking for ramp inputs and step responses

## 🛠️ Technical Implementation

All assignments were implemented in **MATLAB** using:
- Fuzzy Logic Toolbox
- Control System Toolbox
- ANFIS for TSK model training
- Custom MATLAB scripts for simulation and analysis

# computational-intelligence---2025
