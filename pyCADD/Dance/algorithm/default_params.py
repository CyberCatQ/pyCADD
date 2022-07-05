# Default Hyperparameter Grids of Machine Learning Models

GBT_DEFAULT_PARAMS =  {
            # 'gpu_id': [0],
            # 'tree_method': ['gpu_hist'],
            'n_estimators': [200, 250, 300],
            'max_depth': [3, 5, 10, 20],
            'gamma': [0.01, 0.1, 0.5, 1],
            'learning_rate': [0.01, 0.05, 0.1],
            'subsample': [0.3, 0.5, 0.6],
            'alpha': [0.01, 0.1, 0.5, 1],
            'colsample_bytree': [0.3, 0.5, 1],
            'objective': ['binary:logistic'],
            'eval_metric': ['auc'],
            'use_label_encoder': [False]
        }

LR_DEFAULT_PARAMS = {
            'C': [0.01, 0.1, 1],
            'penalty': ['l1', 'l2']
        }

RF_DEFAULT_PARAMS = {
            'n_estimators': [200, 250, 300],
            'max_depth': [3, 5, 10, 20],
            'max_features': ['sqrt', 'log2'],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 5],
        }
