import React, { useEffect } from 'react';
import { initializeApp } from 'firebase/app';
import * as firebaseui from 'firebaseui';
import { getAuth, onAuthStateChanged, signOut } from 'firebase/auth';
import 'firebaseui/dist/firebaseui.css';
import { useNavigate } from 'react-router-dom'; // <-- Import this

// Firebase configuration
const config = {
  apiKey: "AIzaSyB6nkVvXFrHv8lwEzdxj_Q-eWQpmhmtMLE",
  authDomain: "m24convention.firebaseapp.com",
  projectId: "m24convention",
  storageBucket: "m24convention.appspot.com",
  messagingSenderId: "354807424077",
  appId: "1:354807424077:web:fd7ebac2887f21bb17e804",
  measurementId: "G-L17QS3GJZW"
};

const app = initializeApp(config);
const authInstance = getAuth(app);

function startFirebaseUI(containerId, onAuthenticatedCallback, delay = 50) {
  setTimeout(() => {
    let ui = firebaseui.auth.AuthUI.getInstance();
    if (!ui) {
      ui = new firebaseui.auth.AuthUI(authInstance);
    }
    const localUiConfig = {
      ...uiConfig,
      callbacks: {
        signInSuccessWithAuthResult: (authResult, redirectUrl) => {
          if (onAuthenticatedCallback) {
            onAuthenticatedCallback(authResult.user);
          }
          return false;
        },
      },
    };
    ui.start(containerId, localUiConfig);
  }, delay);
}

function deleteFirebaseUI() {
  const existingUi = firebaseui.auth.AuthUI.getInstance();
  if (existingUi) {
    existingUi.delete();
  }
}

const onAuthChange = (callback) => {
  return onAuthStateChanged(authInstance, callback);
};

const uiConfig = {
  signInOptions: [
    'google.com',
    'password',
  ],
};
function Login() {
  const navigate = useNavigate(); // <-- Use the hook here
  
  useEffect(() => {
    startFirebaseUI('#firebaseui-auth-container', (user) => {
      navigate('/');  // <-- Navigate after successful login
    });

    return () => {
      deleteFirebaseUI();
    };
  }, [navigate]);

  return <div id="firebaseui-auth-container"></div>;
}

export { Login, authInstance, onAuthChange, signOut };
