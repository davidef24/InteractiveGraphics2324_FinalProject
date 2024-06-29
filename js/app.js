// Constants
import * as THREE from './three.module.js'

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000 );


const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setAnimationLoop(animate);
document.body.appendChild(renderer.domElement);

const geometry = new THREE.BoxGeometry(7, 4, 4);

// Adjust material properties for transparency and wireframe
const material = new THREE.MeshBasicMaterial({
    color: 0xFAEBD7,
    transparent: true,
    wireframe: true,
    wireframeLinewidth: 2,
});

const cube = new THREE.Mesh(geometry, material);
scene.add(cube);

camera.position.z = 5;

let isRotating = false;
let previousMousePosition = {
    x: 0,
    y: 0
};

// Function to handle mouse down event
function handleMouseDown(event) {
    isRotating = true;
    previousMousePosition = {
        x: event.clientX,
        y: event.clientY
    };
}

function handleResize(event) {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

// Function to handle mouse move event
function handleMouseMove(event) {
    if (isRotating) {
        let deltaMove = {
            x: event.clientX - previousMousePosition.x,
            y: event.clientY - previousMousePosition.y
        };

        cube.rotation.y += deltaMove.x * 0.01;
        cube.rotation.x += deltaMove.y * 0.01;

        previousMousePosition = {
            x: event.clientX,
            y: event.clientY
        };
    }
}

// Function to handle mouse wheel event (zooming)
function handleMouseWheel(event) {
    let delta = event.deltaY * 0.001;
    camera.position.z += delta;
}

// Function to handle mouse up event
function handleMouseUp(event) {
    isRotating = false;
}

// Add mouse event listeners
renderer.domElement.addEventListener('mousedown', handleMouseDown, false);
renderer.domElement.addEventListener('mousemove', handleMouseMove, false);
renderer.domElement.addEventListener('mouseup', handleMouseUp, false);
renderer.domElement.addEventListener('wheel', handleMouseWheel, false);
window.addEventListener('resize', handleResize);


function animate() {
    renderer.render(scene, camera);
}
