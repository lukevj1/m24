import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import { Calendar, momentLocalizer, Views } from 'react-big-calendar';
import { getFirestore, collection, onSnapshot, query } from 'firebase/firestore';
import moment from 'moment';
import { useSwipeable } from 'react-swipeable';
import './schedule.css';

import 'react-big-calendar/lib/css/react-big-calendar.css';

const localizer = momentLocalizer(moment);

function Schedule() {
    const [events, setEvents] = useState([]);
    const [currentDate, setCurrentDate] = useState(new Date(2024, 6, 26));
    const navigate = useNavigate();
    useEffect(() => {
        const db = getFirestore();
        const eventsRef = collection(db, 'events');

        const q = query(eventsRef);

        const unsubscribe = onSnapshot(q, snapshot => {
            const fetchedEvents = snapshot.docs.map(doc => {
                const data = doc.data();
                return {
                    id: doc.id,
                    title: data.title,
                    start: new Date(data.start_time.seconds * 1000),
                    end: new Date((data.start_time.seconds + data.duration * 3600) * 1000),
                    desc: data.description,
                    location: data.location,
                    fileURL: data.fileURL
                };
            });
            setEvents(fetchedEvents);
        });

        return () => unsubscribe();
    }, []);

    const navigateDate = (increment) => {
        const newDate = new Date(currentDate);
        newDate.setDate(currentDate.getDate() + increment);
        if (newDate >= new Date(2024, 6, 26) && newDate <= new Date(2024, 6, 30)) {
            setCurrentDate(newDate);
        }
    };
    const handleEventClick = (event) => {
        navigate(`/event/${event.id}`);
    };
    
    const handlers = useSwipeable({
        onSwipedLeft: () => navigateDate(1),
        onSwipedRight: () => navigateDate(-1),
        preventDefaultTouchmoveEvent: true,
        trackMouse: true
    });

    const eventsOnCurrentDay = events.filter(event => moment(event.start).isSame(currentDate, 'day'));
    
    let minTime = new Date(currentDate);
    let maxTime = new Date(currentDate);
    if (eventsOnCurrentDay.length > 0) {
        const firstEventTime = Math.min(...eventsOnCurrentDay.map(event => event.start));
        const lastEventEndTime = Math.max(...eventsOnCurrentDay.map(event => event.end));
        minTime.setHours(new Date(firstEventTime).getHours() - 0.5);
        maxTime.setHours(new Date(lastEventEndTime).getHours() + 0.5);

        if (maxTime.getHours() >= 23.5) {
            maxTime.setHours(24);
        }
    }
    const customEventStyle = (event, start, end, isSelected) => {
        const backgroundImageStyle = event.fileURL 
            ? { 
                backgroundImage: `url(${event.fileURL})`, 
                backgroundSize: 'cover', 
                backgroundPosition: 'center',
                color: 'white',
                // A generic inset shadow to simulate image edge extension:
                boxShadow: 'inset 0 0 50px rgba(0, 0, 0, 0.5)'
              } 
            : {};
    
        return {
            style: {
                ...backgroundImageStyle
            }
        };
    };
    
    const totalDays = 5;
    const currentDay = currentDate.getDate() - 25; // Subtracting the start day to get the current index
    return (
        <div>
            <h2>Conference Schedule</h2>
            <h3>{moment(currentDate).format('MMMM D, YYYY')}</h3>
            <div style={{ display: 'flex', alignItems: 'center' }} {...handlers}>
                <Calendar
                    localizer={localizer}
                    events={events}
                    views={[Views.DAY]}
                    date={currentDate}
                    onNavigate={date => setCurrentDate(date)}
                    defaultView={Views.DAY}
                    startAccessor="start"
                    endAccessor="end"
                    timeslots={1}
                    min={minTime}
                    max={maxTime}
                    style={{ height: "500px", flex: 1 }}
                    toolbar={false}
                    eventPropGetter={customEventStyle}
                    onSelectEvent={handleEventClick}
                />
            </div>
            {eventsOnCurrentDay.length === 0 && <p>No events scheduled for this day.</p>}
            <div style={{ display: 'flex', justifyContent: 'center', marginTop: '10px' }}>
                {[...Array(totalDays).keys()].map(day => (
                    <div 
                        key={day}
                        style={{
                            width: '20px', 
                            height: '5px', 
                            margin: '0 3px', 
                            backgroundColor: day < currentDay ? 'blue' : 'lightgray'
                        }}
                    />
                ))}
            </div>
        </div>
    );
}

export default Schedule;
