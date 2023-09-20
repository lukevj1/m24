import React, { useEffect, useState } from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import { Login, authInstance, onAuthChange, signOut } from './components/firebase';
import Schedule from './components/schedule';
import EventDetail from './components/eventDetail'; 
import './App.css';

function PrivateRoute({ children, user }) {  // Added the 'user' prop
  if (user) {
    return children; 
  } else {
    return <Navigate to="/login" replace />;
  }
}

function App() {
  const [user, setUser] = useState(null);

  useEffect(() => {
    const unsubscribe = onAuthChange(setUser);
    return () => unsubscribe();
  }, []);

  const handleSignOut = async () => {
    try {
      await signOut(authInstance);
      setUser(null);
    } catch (error) {
      console.error("Error signing out:", error);
    }
  };


  return (
    <Router>
      <div id="app">
        <header style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '10px 20px', boxShadow: '0px 2px 10px rgba(0,0,0,0.1)' }}>
          <div style={{ display: 'flex', alignItems: 'center' }}>
            <img src={process.env.PUBLIC_URL + '/logo.png'} alt="M24 Convention Logo" id="logo" style={{ marginRight: '20px' }} />
            <h1>M24 Convention</h1>
          </div>
          {user ? (
            <div style={{ display: 'flex', alignItems: 'center' }}>
              <p style={{ margin: '0 10px' }}>Welcome, {user.displayName || user.email}</p>
              <button onClick={handleSignOut} style={{ padding: '5px 15px', borderRadius: '5px', border: 'none', background: '#333', color: '#fff' }}>Sign Out</button>
            </div>
          ) : (
            <p>The future is now</p>
          )}
        </header>
        <section>
        <Routes>
          <Route path="/login" element={<Login />} />
          <Route path="/" element={
            <PrivateRoute user={user}>
              <Schedule />
            </PrivateRoute>
          }/>
          <Route path="/event/:eventId" element={
            <PrivateRoute user={user}>
              <EventDetail />
            </PrivateRoute>
          }/>
        </Routes>
        </section>
      </div>
    </Router>
  );
}

export default App;
