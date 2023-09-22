import React from 'react';
import { BrowserRouter as Router, Routes, Route, Navigate } from 'react-router-dom';
import Schedule from './components/schedule';
import EventDetail from './components/eventDetail';
import './App.css';
import { app } from './components/firebase';

function App() {
  return (
    <Router>
      <div id="app">
        <header style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '10px 20px', boxShadow: '0px 2px 10px rgba(0,0,0,0.1)' }}>
          <div style={{ display: 'flex', alignItems: 'center' }}>
            <img src={process.env.PUBLIC_URL + '/logo.png'} alt="M24 Convention Logo" id="logo" style={{ marginRight: '20px' }} />
            <h1>M24 Convention</h1>
          </div>
          <p>The future is now</p>
        </header>
        <section>
          <Routes>
            <Route path="/" element={<Navigate to="/schedule" replace />} />
            <Route path="/schedule" element={<Schedule />} />
            <Route path="/event/:eventId" element={<EventDetail />} />
          </Routes>
        </section>
      </div>
    </Router>
  );
}

export default App;
