#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick Start Example - Advanced Dual-Target Affinity Prediction

This script demonstrates the basic usage of the system.
Run this to verify installation and see the system in action.
"""

import sys
import numpy as np

def check_dependencies():
    """Check if all required packages are installed"""
    print("Checking dependencies...")
    
    required_packages = {
        'torch': 'PyTorch',
        'numpy': 'NumPy',
        'pandas': 'Pandas',
        'sklearn': 'scikit-learn',
        'scipy': 'SciPy',
        'matplotlib': 'Matplotlib',
        'seaborn': 'Seaborn',
        'tqdm': 'tqdm'
    }
    
    missing = []
    for package, name in required_packages.items():
        try:
            __import__(package)
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  ✗ {name} (MISSING)")
            missing.append(name)
    
    if missing:
        print(f"\n❌ Missing packages: {', '.join(missing)}")
        print("\nInstall with:")
        print("  pip install torch numpy pandas scikit-learn scipy matplotlib seaborn tqdm")
        return False
    else:
        print("\n✅ All dependencies installed!")
        return True


def quick_demo():
    """
    Quick demonstration of the system
    
    This creates a minimal example showing the workflow
    """
    print("\n" + "="*70)
    print("Quick Demo - Advanced Dual-Target Affinity Prediction")
    print("="*70)
    
    try:
        # Import modules
        print("\n1. Importing modules...")
        from advanced_dual_target_mtl import (
            AdvancedConfig,
            AdvancedDualTargetPredictor,
            AdvancedTrainer,
            generate_synthetic_data,
            prepare_data_loaders,
            evaluate_model,
            set_seed
        )
        print("   ✓ Modules imported successfully")
        
        # Configuration
        print("\n2. Setting up configuration...")
        config = AdvancedConfig()
        config.epochs = 10  # Reduced for demo
        config.batch_size = 32
        config.cv_folds = 2  # Reduced for demo
        set_seed(config.random_seed)
        print("   ✓ Configuration ready")
        
        # Generate data
        print("\n3. Generating synthetic data...")
        drug_features, protein1_features, protein2_features, \
            labels_task1, labels_task2 = generate_synthetic_data(
                n_samples=200,  # Small dataset for demo
                drug_dim=config.drug_feature_dim,
                protein_dim=config.protein_feature_dim,
                random_seed=config.random_seed
            )
        print(f"   ✓ Generated {len(drug_features)} samples")
        
        # Prepare data loaders
        print("\n4. Preparing data loaders...")
        train_loader, val_loader, test_loader, scalers = prepare_data_loaders(
            drug_features, protein1_features, protein2_features,
            labels_task1, labels_task2, config,
            train_ratio=0.7, val_ratio=0.15
        )
        print("   ✓ Data loaders ready")
        
        # Create model
        print("\n5. Creating model...")
        model = AdvancedDualTargetPredictor(config)
        print(f"   ✓ Model created with {model.count_parameters():,} parameters")
        
        # Create trainer
        print("\n6. Creating trainer...")
        trainer = AdvancedTrainer(model, config)
        print("   ✓ Trainer ready")
        
        # Train (just a few epochs for demo)
        print("\n7. Training model (demo: 10 epochs)...")
        print("   This may take a few minutes...")
        history = trainer.train(train_loader, val_loader)
        print("   ✓ Training completed")
        
        # Evaluate
        print("\n8. Evaluating model...")
        test_metrics = evaluate_model(model, test_loader, config.device)
        print("   ✓ Evaluation completed")
        
        # Results
        print("\n" + "="*70)
        print("Demo Results")
        print("="*70)
        print(f"\nTask 1:")
        print(f"  RMSE: {test_metrics['task1_rmse']:.4f}")
        print(f"  R²:   {test_metrics['task1_r2']:.4f}")
        print(f"\nTask 2:")
        print(f"  RMSE: {test_metrics['task2_rmse']:.4f}")
        print(f"  R²:   {test_metrics['task2_r2']:.4f}")
        
        print("\n" + "="*70)
        print("✅ Demo completed successfully!")
        print("="*70)
        print("\nNext steps:")
        print("  1. Run full training: python advanced_main.py")
        print("  2. Check outputs/ folder for results")
        print("  3. Check figures/ folder for visualizations")
        print("  4. Read README.md for more information")
        print("="*70 + "\n")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Error during demo: {e}")
        print("\nTroubleshooting:")
        print("  1. Make sure all dependencies are installed")
        print("  2. Check training.log for detailed error messages")
        print("  3. See INSTALLATION_GUIDE.md for help")
        return False


def main():
    """Main function"""
    print("\n" + "="*70)
    print("Advanced Dual-Target Affinity Prediction - Quick Start")
    print("="*70)
    
    # Check dependencies
    if not check_dependencies():
        return 1
    
    # Run demo
    print("\nDo you want to run a quick demo? (y/n): ", end='')
    try:
        response = input().strip().lower()
    except (EOFError, KeyboardInterrupt):
        print("\n\nDemo skipped.")
        return 0
    
    if response in ['y', 'yes']:
        success = quick_demo()
        return 0 if success else 1
    else:
        print("\nDemo skipped.")
        print("\nTo run full training:")
        print("  python advanced_main.py")
        print("\nFor more information:")
        print("  See README.md")
        return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\nInterrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nUnexpected error: {e}")
        sys.exit(1)
