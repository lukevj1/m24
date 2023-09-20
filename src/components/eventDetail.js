import React, { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { getFirestore, doc, getDoc } from 'firebase/firestore';
import { useNavigate } from 'react-router-dom';

function EventDetail() {
    const { eventId } = useParams();
    const [eventData, setEventData] = useState(null);
    const navigate = useNavigate();

    useEffect(() => {
        const fetchData = async () => {
            const db = getFirestore();
            const eventRef = doc(db, 'events', eventId);
            const eventSnapshot = await getDoc(eventRef);
            
            if (eventSnapshot.exists()) {
                setEventData(eventSnapshot.data());
            }
        };
        
        fetchData();
    }, [eventId]);

    if (!eventData) return <div>Loading...</div>;

    return (
        <div id="app">
            <button style={{marginBottom: '1rem'}} onClick={() => navigate(-1)}>Back</button>
            <section>
                <header>
                    <h1>{eventData.title}</h1>
                </header>
                {/* Display other event details as desired in a modern layout */}
                {/* Example: */}
                <img id="logo" src={eventData.fileURL} alt="Event" />
                <p>{eventData.description}</p>
                <p>Location: {eventData.location}</p>
                {/* Add more details as needed */}
            </section>
        </div>
    );
}

export default EventDetail;