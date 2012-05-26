/** demos.h
	Copyright (C) 2012 Jackson Lee

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/

#include "btBulletDynamicsCommon.h"
#include "Fluids/fluid_system.h"
#include "Fluids/SPHInterface.h"


class FluidSystemDemo
{
protected:
	btAlignedObjectArray<btRigidBody*> m_rigidBodies;

public:
	virtual ~FluidSystemDemo() {}
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes) {}
	virtual void deactivate() 
	{
		for(int i = 0; i < m_rigidBodies.size(); ++i)
		{
			if(m_rigidBodies[i])
			{
				if( m_rigidBodies[i]->getMotionState() ) delete m_rigidBodies[i]->getMotionState();
				delete m_rigidBodies[i];
			}
		}
		
		m_rigidBodies.clear();
	}
	
	virtual void addToWorld(btDynamicsWorld *world) { for(int i = 0; i < m_rigidBodies.size(); ++i)world->addRigidBody(m_rigidBodies[i]); }
	virtual void removeFromWorld(btDynamicsWorld *world) { for(int i = 0; i < m_rigidBodies.size(); ++i)world->removeRigidBody(m_rigidBodies[i]); }
	
	//	rename as 'stepSimulation'?
	virtual void update(FluidSystem *fluidSystem) {}
	
	virtual void reset(FluidSystem *fluidSystem, int maxFluidParticles) = 0;
};
class Demo_DamBreak : public FluidSystemDemo
{
public:
	virtual void reset(FluidSystem *fluidSystem, int maxFluidParticles)
	{
		const float VOL_BOUND = 20.0f;
		Vector3DF volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		Vector3DF volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluidSystem->initialize(maxFluidParticles, volMin, volMax);
		
		const float INIT_BOUND = 20.0f;
		Vector3DF initMin(-INIT_BOUND, 10.f, 0.f);
		Vector3DF initMax(INIT_BOUND, 55.f, INIT_BOUND);
		FluidEmitter::addVolume( fluidSystem, initMin, initMax, fluidSystem->getEmitterSpacing() * 0.87 );
	}
};
class Demo_Drop : public FluidSystemDemo
{
public:
	virtual void reset(FluidSystem *fluidSystem, int maxFluidParticles)
	{
		const float VOL_BOUND = 30.0f;
		Vector3DF volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		Vector3DF volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluidSystem->initialize(maxFluidParticles, volMin, volMax);
	
		const float INIT_BOUND = 20.0f;
		Vector3DF initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		Vector3DF initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluidSystem, initMin, initMax, fluidSystem->getEmitterSpacing() * 0.87 );
	}
};
class Demo_EmitterAndAbsorber : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 0.0;
		const btVector3 LOCAL_INERTIA(0,0,0);	//Calculated only if mass != 0

		btCollisionShape *wallShape = new btBoxShape( btVector3(5.0, 15.0, 50.0) );
		collisionShapes->push_back(wallShape);
		
		btTransform transform;
		btDefaultMotionState *motionState;
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(25.0, 10.0, 0.0) );
		motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfoA(MASS, motionState, wallShape, LOCAL_INERTIA);
		m_rigidBodies.push_back( new btRigidBody(rbInfoA) );
		
		transform = btTransform( btQuaternion::getIdentity(), btVector3(-25.0, 10.0, 0.0) );
		motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfoB(MASS, motionState, wallShape, LOCAL_INERTIA);
		m_rigidBodies.push_back( new btRigidBody(rbInfoB) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(0.0, 10.0, 25.0) );
		motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfoC(MASS, motionState, wallShape, LOCAL_INERTIA);
		m_rigidBodies.push_back( new btRigidBody(rbInfoC) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(55.0, 10.0, -25.0) );
		motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfoD(MASS, motionState, wallShape, LOCAL_INERTIA);
		m_rigidBodies.push_back( new btRigidBody(rbInfoD) );
		
		transform = btTransform( btQuaternion(SIMD_PI*0.5, 0.0, 0.0), btVector3(-55.0, 10.0, -25.0) );
		motionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfoE(MASS, motionState, wallShape, LOCAL_INERTIA);
		m_rigidBodies.push_back( new btRigidBody(rbInfoE) );
	}
	
	virtual void update(FluidSystem *fluidSystem)
	{
		FluidEmitter emitter;
		emitter.m_pitch = -90.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(10.f, 10.f, 10.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 20;
			emitter.emit( fluidSystem, NUM_PARTICLES_EMITTED, fluidSystem->getEmitterSpacing() );
		}
		
		//
		FluidAbsorber absorber;
		absorber.m_min.setValue(-100.0, -20.0, -75.0);
		absorber.m_max.setValue(100.0, 20.0, -50.0);
		absorber.absorb(fluidSystem);
	}
	
	virtual void reset(FluidSystem *fluidSystem, int maxFluidParticles)
	{
		const float VOL_BOUND = 120.0f;
		Vector3DF volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		Vector3DF volMax(VOL_BOUND, VOL_BOUND, VOL_BOUND);
		
		fluidSystem->initialize(maxFluidParticles, volMin, volMax);
	
		//const float INIT_BOUND = 20.0f;
		//Vector3DF initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		//Vector3DF initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		//FluidEmitter::addVolume( fluidSystem, initMin, initMax, fluidSystem->getEmitterSpacing() * 0.87 );
	}
};
class Demo_DynamicBox : public FluidSystemDemo
{
public:
	virtual void initialize(btAlignedObjectArray<btCollisionShape*> *collisionShapes)
	{
		const btScalar MASS = 1.f;
	
		btCollisionShape* largeBoxShape = new btBoxShape( btVector3(10.0, 10.0, 10.0) );
		collisionShapes->push_back(largeBoxShape);

		//Rigid bodies are dynamic if and only if mass is non zero, otherwise static
		btVector3 localInertia(0,0,0);
		if(MASS != 0.f) largeBoxShape->calculateLocalInertia(MASS, localInertia);
		
		btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 10.0, 0.0) );

		btDefaultMotionState* motionState = new btDefaultMotionState(startTransform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(MASS, motionState, largeBoxShape, localInertia);
		m_rigidBodies.push_back( new btRigidBody(rbInfo) );
	}
	
	virtual void update(FluidSystem *fluidSystem)
	{
		FluidEmitter emitter;
		emitter.m_pitch = -90.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(10.f, 10.f, 10.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 5;
			emitter.emit( fluidSystem, NUM_PARTICLES_EMITTED, fluidSystem->getEmitterSpacing() );
		}
	}
	
	virtual void reset(FluidSystem *fluidSystem, int maxFluidParticles)
	{
		const float VOL_BOUND = 30.0f;
		Vector3DF volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		Vector3DF volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		fluidSystem->initialize(maxFluidParticles, volMin, volMax);
		
		const float INIT_BOUND = 20.0f;
		Vector3DF initMin(-INIT_BOUND, 20.0f, -INIT_BOUND);
		Vector3DF initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( fluidSystem, initMin, initMax, fluidSystem->getEmitterSpacing() * 0.87 );
	
		if( m_rigidBodies.size() && m_rigidBodies[0] )
		{
			btTransform startTransform( btQuaternion::getIdentity(), btVector3(0.0, 10.0, 0.0) );
			m_rigidBodies[0]->setCenterOfMassTransform(startTransform);
		}
	}
};

	//	next: add a bucket composed of constrained btBoxShapes
class Demo_Bucket : public FluidSystemDemo
{
public:
};
	
