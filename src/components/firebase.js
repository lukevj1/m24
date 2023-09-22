import { initializeApp } from 'firebase/app';

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

// Export only necessary items
export { app };
